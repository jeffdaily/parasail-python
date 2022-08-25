import codecs
import os
import platform
import re
import shutil
import stat
import subprocess
import time
try:
    from urllib import urlretrieve
except ImportError:
    from urllib.request import urlretrieve
import tarfile
import zipfile

from distutils.dir_util import copy_tree
from distutils.util import get_platform
from setuptools import setup
from setuptools.command.install import install as install_

###############################################################################

NAME = "parasail"
PACKAGES = ['parasail']
META_PATH = os.path.join("parasail", "__init__.py")
KEYWORDS = ["Smith-Waterman", "Needleman-Wunsch"]
CLASSIFIERS = [
    "Development Status :: 5 - Production/Stable",
    "Intended Audience :: Science/Research",
    "Natural Language :: English",
    "License :: OSI Approved :: BSD License",
    "Operating System :: OS Independent",
    "Programming Language :: C",
    "Programming Language :: Python",
    "Programming Language :: Python :: 2",
    "Programming Language :: Python :: 2.5",
    "Programming Language :: Python :: 2.6",
    "Programming Language :: Python :: 2.7",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.0",
    "Programming Language :: Python :: 3.1",
    "Programming Language :: Python :: 3.2",
    "Programming Language :: Python :: 3.3",
    "Programming Language :: Python :: 3.4",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Topic :: Software Development :: Libraries :: Python Modules",
]
INSTALL_REQUIRES = ["numpy"]

###############################################################################

TIMEOUT = 10
HERE = os.path.abspath(os.path.dirname(__file__))


def is_python_64bit():
    import sys
    return sys.maxsize > 2**32

def read(*parts):
    """
    Build an absolute path from *parts* and and return the contents of the
    resulting file.  Assume UTF-8 encoding.
    """
    with codecs.open(os.path.join(HERE, *parts), "rb", "utf-8") as f:
        return f.read()


META_FILE = read(META_PATH)


def find_meta(meta):
    """
    Extract __*meta*__ from META_FILE.
    """
    meta_match = re.search(
        r"^__{meta}__ = ['\"]([^'\"]*)['\"]".format(meta=meta),
        META_FILE, re.M
    )
    if meta_match:
        return meta_match.group(1)
    raise RuntimeError("Unable to find __{meta}__ string.".format(meta=meta))


URI = find_meta("uri")

def get_libname():
    libname = "libparasail.so"
    if platform.system() == "Darwin":
        libname = "libparasail.dylib"
    elif platform.system() == "Windows":
        libname = "parasail.dll"
    return libname

def get_includes():
    return ["include/*.h",
            "include/parasail/*.h",
            "include/parasail/matrices/*.h"]

def unzip(archive, destdir):
    thefile=zipfile.ZipFile(archive)
    thefile.extractall(destdir)
    thefile.close()

def untar(archive):
    tar = tarfile.open(archive)
    tar.extractall()
    tar.close()

def find_file(filename, start="."):
    for root, dirs, files in os.walk(start, topdown=False):
        for name in files:
            if name == filename:
                return root
    return None

# attempt to run parallel make with at most 8 workers
def cpu_count():
    try:
        import multiprocessing
        return min(8, multiprocessing.cpu_count())
    except:
        return 1

# unzipping parasail C library zip file does not preserve executable permissions
def fix_permissions(start):
    execmode = stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
    filenames = [
        "version.sh",
        "func_group_rowcols.py",
        "codegen.py",
        "satcheck.py",
        "func_group_tables.py",
        "tester.py",
        "diff_all.sh",
        "names.py",
        "dispatcher.py",
        "isastubs.py",
        "funcs.py",
        "makedef.py",
        "gap_tester.py",
        "func_groups.py",
        "configure",
        "install-sh",
        "config.sub",
        "test-driver",
        "config.guess",
        "compile",
        "missing",
        "depcomp",
        ]
    for root, dirs, files in os.walk(start, topdown=False):
        for name in files:
            if name in filenames:
                fullpath = os.path.join(root,name)
                st = os.stat(fullpath)
                os.chmod(fullpath, st.st_mode | execmode)

def run_autoreconf(root):
    print("Running autoreconf -fi from {}".format(root))
    all_good = True
    tools = ['perl', 'm4', 'autoconf', 'automake', 'libtoolize', 'autoreconf']
    if platform.system() == "Darwin":
        tools[4] = 'glibtoolize'
        os.environ['LIBTOOL'] = 'glibtool'
        os.environ['LIBTOOLIZE'] = 'glibtoolize'
    for tool in tools:
        try:
            output = subprocess.check_output([tool, '--version'],
                    stderr=subprocess.STDOUT)
            print(output.strip().splitlines()[0])
        except subprocess.CalledProcessError as e:
            print('{} --version failed'.format(tool))
            all_good = False
        except OSError as e:
            print('{} not found'.format(tool))
            all_good = False
    if all_good:
        retcode = -1
        try:
            retcode = subprocess.Popen([
                'autoreconf',
                '-fi',
                ], cwd=root).wait()
        except Exception as e:
            all_good = False
        else:
            print("autoreconf -fi exited with return code {}".format(retcode))
            if retcode != 0:
                all_good = False
    return all_good

def build_autotools(patch_m4=False):
    if patch_m4:
        print("Building autotools with patched m4")
    else:
        print("Building autotools")
    save_cwd = os.getcwd()
    top = os.path.join(os.getcwd(), 'autotools')
    if not os.path.exists(top):
        os.mkdir(top)
    # we know these versions to work
    tools = [('m4', '1.4.17'),
            ('autoconf', '2.69'),
            ('automake', '1.15.1'),
            ('libtool', '2.4.6')]
    for tool,version in tools:
        os.chdir(top)
        tdir = '{}-{}'.format(tool, version)
        tarball = '{}.tar.gz'.format(tdir)
        binary = '{}/bin/{}'.format(top, tool)
        url = 'http://ftp.gnu.org/gnu/{}/{}'.format(tool, tarball)
        if os.path.exists(tarball):
            print("{} already exists! Using existing copy.".format(tarball))
        else:
            print("Downloading {}".format(url))
            for attempt in range(10):
                try:
                    name,hdrs = urlretrieve(url, tarball)
                except Exception as e:
                    print(repr(e))
                    print("Will retry in {} seconds".format(TIMEOUT))
                    time.sleep(TIMEOUT)
                else:
                    break
            else:
                # we failed all the attempts - deal with the consequences.
                print("All attempts to download {} have failed".format(tarball))
                os.chdir(save_cwd)
                return False
        if os.path.exists(tdir):
            print("{} already exists! Using existing sources.".format(tdir))
        else:
            print("Expanding {}".format(tarball))
            untar(tarball)
        if os.path.exists(binary):
            print("{} already exists! Skipping build.".format(binary))
        else:
            os.chdir(os.path.join(top,tdir))
            CPPFLAGS=""
            if tool == 'm4' and patch_m4:
                CPPFLAGS="CPPFLAGS=-D_IO_IN_BACKUP=0x100"
                for filename in ['lib/freading.c', 'lib/fseeko.c', 'lib/fpurge.c', 'lib/freadahead.c', 'lib/fflush.c']:
                    with open(filename, "r") as source:
                        lines = source.readlines()
                    newlines = [line.replace('_IO_ftrylockfile','_IO_EOF_SEEN') for line in lines]
                    if lines == newlines:
                        print("{} skipped".format(filename))
                    else:
                        with open(filename, "w") as source:
                            for line in newlines:
                                source.write(line)
                        print("{} patched".format(filename))
            print("configuring {}".format(tool))
            print("running ./configure --prefix={} {}".format(top,CPPFLAGS))
            proc = subprocess.Popen([
                './configure', '--prefix={}'.format(top), '{}'.format(CPPFLAGS)
                ], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            stdout,stderr = proc.communicate()
            if 0 != proc.returncode:
                print(stdout)
                print("configure of {} failed".format(tool))
                os.chdir(save_cwd)
                return False
            print("making and installing {}".format(tool))
            proc = subprocess.Popen([
                'make', '-j', str(cpu_count()), 'install'
                ], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            stdout,stderr = proc.communicate()
            if 0 != proc.returncode:
                print(stdout)
                print("make of {} failed".format(tool))
                os.chdir(save_cwd)
                return False
    os.chdir(save_cwd)
    return True

# Download, unzip, configure, and make parasail C library from github.
# Attempt to skip steps that may have already completed.
def build_parasail(libname):
    archive = 'parasail-master.zip'
    unzipped_archive = 'parasail-master'
    destdir = os.getcwd()

    if not os.path.exists(archive):
        print("Downloading latest parasail master")
        theurl = 'https://github.com/jeffdaily/parasail/archive/master.zip'
        for attempt in range(10):
            try:
                name,hdrs = urlretrieve(theurl, archive)
            except Exception as e:
                print(repr(e))
                print("Will retry in {} seconds".format(TIMEOUT))
                time.sleep(TIMEOUT)
            else:
                break
        else:
            # we failed all the attempts - deal with the consequences.
            raise RuntimeError("All attempts to download latest parasail master have failed")
    else:
        print("Archive '{}' already downloaded".format(archive))

    if not os.path.exists(unzipped_archive):
        print("Unzipping parasail master archive")
        unzip(archive, destdir)
    else:
        print("Archive '{}' already unzipped to {}".format(archive,destdir))

    # need to search for a file specific to top-level parasail archive
    parasail_root = find_file('version.sh')
    if not os.access(os.path.join(parasail_root,'version.sh'), os.X_OK):
        print("fixing executable bits after unzipping")
        fix_permissions(parasail_root)
    else:
        print("parasail archive executable permissions ok")

    root = find_file('configure', parasail_root)
    if root is None:
        print("Unable to find parasail configure script")
        root = find_file('configure.ac', parasail_root)
        # force use of /usr/bin/perl
        # on OSX, homebrew perl might not support threading and automake
        # build will fail
        print('Prepending /usr/bin to PATH')
        os.environ['PATH'] = '/usr/bin' + os.pathsep + os.environ['PATH']
        print("PATH={}".format(os.environ['PATH']))
        if not run_autoreconf(root):
            newpath = os.path.join(os.getcwd(), 'autotools', 'bin')
            print("Prepending {} to PATH".format(newpath))
            os.environ['PATH'] = newpath + os.pathsep + os.environ['PATH']
            if platform.system() == "Darwin" and 'M4' not in os.environ:
                os.environ['M4'] = '/usr/bin/m4'
            print("PATH={}".format(os.environ['PATH']))
            if not build_autotools():
                if not build_autotools(True):
                    raise RuntimeError("building autotools failed")
            if not run_autoreconf(root):
                raise RuntimeError("autoreconf -fi failed")
    root = find_file('configure', parasail_root)
    if root is None:
        raise RuntimeError("Unable to find parasail configure script after autoreconf step")

    if find_file('config.status', parasail_root) is None:
        print("configuring parasail in directory {}".format(parasail_root))
        retcode = subprocess.Popen([
            './configure',
            '--enable-shared',
            '--disable-static'
            ], cwd=parasail_root).wait()
        if 0 != retcode:
            raise RuntimeError("configure failed")
    else:
        print("parasail already configured in directory {}".format(parasail_root))

    if find_file(libname, parasail_root) is None:
        print("making parasail")
        retcode = subprocess.Popen([
            'make',
            '-j',
            str(cpu_count())
            ], cwd=parasail_root).wait()
        if 0 != retcode:
            raise RuntimeError("make failed")
    else:
        print("parasail library '{}' already made".format(libname))
    src = os.path.join(parasail_root, '.libs', libname)
    dst = 'parasail'
    print("copying {} to {}".format(src,dst))
    shutil.copy(src,dst)

    # copy headers into the parasail directory, so they can be added to the wheel
    parasail_headers = os.path.join(parasail_root, "parasail")
    dst = os.path.join("parasail", "include", "parasail")
    print("copying {} to {}".format(parasail_headers, dst))
    copy_tree(parasail_headers, dst)
    parasail_h = os.path.join(parasail_root, "parasail.h")
    dst = os.path.join("parasail", "include")
    print("copying {} to {}".format(parasail_h, dst))
    shutil.copy(parasail_h, dst)

def github_api_json(address):
    import json
    import sys
    if (sys.version_info > (3, 0)):
        import urllib.request
        with urllib.request.urlopen(address) as url:
            data = json.loads(url.read().decode())
    else:
        import urllib
        response = urllib.urlopen(address)
        data = json.loads(response.read())
    return data

def download_windows_dll():
    libname = get_libname()
    libpath = os.path.join("parasail", libname)
    print("Downloading latest parasail release info from github")
    address = "https://api.github.com/repos/jeffdaily/parasail/releases/latest"
    data = None
    for attempt in range(10):
        try:
            data = github_api_json(address)
            if not data or 'assets' not in data:
                raise RuntimeError("Unable to download github asset JSON from "+address)
        except Exception as e:
            print(repr(e))
            print("Will retry in {} seconds".format(TIMEOUT))
            time.sleep(TIMEOUT)
        else:
            break
    else:
        # we failed all the attempts - deal with the consequences.
        raise RuntimeError("All attempts to download github asset JSON have failed")
    asset = None
    search = "win32-v140"
    if is_python_64bit():
        search = "win64-v140"
    for maybe_asset in data['assets']:
        if search in maybe_asset['browser_download_url']:
            asset = maybe_asset['browser_download_url']
            break
    if not asset:
        raise RuntimeError("Unable to determine asset URL")
    print("Downloading latest parasail release {}".format(asset))
    archive = asset.rsplit('/',1)[-1]
    for attempt in range(10):
        try:
            name,hdrs = urlretrieve(asset, archive)
        except Exception as e:
            print(repr(e))
            print("Will retry in {} seconds".format(TIMEOUT))
            time.sleep(TIMEOUT)
        else:
            break
    else:
        # we failed all the attempts - deal with the consequences.
        raise RuntimeError("All attempts to download asset URL have failed")
    unzipped_archive_destination = archive.rsplit('.',1)[0]
    print("Unzipping {}".format(archive))
    unzip(archive, unzipped_archive_destination)
    print("Locating {}".format(libname))
    root = find_file(libname)
    src = os.path.join(root, libname)
    dst = 'parasail'
    print("copying {} to {}".format(src,dst))
    shutil.copy(src,dst)

    # copy headers into the parasail directory, so they can be added to the wheel
    headers_src = os.path.join(root, "..", "include")
    dst = os.path.join("parasail", "include")
    print("copying {} to {}".format(headers_src, dst))
    copy_tree(headers_src, dst)

    # clean up downloaded asset
    shutil.rmtree(unzipped_archive_destination)

def prepare_shared_lib():
    libname = get_libname()
    libpath = os.path.join("parasail", libname)
    if not os.path.exists(libpath):
        if "PARASAIL_SKIP_BUILD" in os.environ:
            print("{} not found, but PARASAIL_SKIP_BUILD set by env var, skipping download and build".format(libpath))
            return
        print("{} not found, attempting to build".format(libpath))
        if platform.system() == "Windows":
            download_windows_dll()
        else:
            build_parasail(libname)
    else:
        print("{} found, skipping build".format(libpath))
    if not os.path.exists(libpath):
        raise RuntimeError("Unable to find shared library {}.".format(libname))

class install(install_):
    def run(self):
        prepare_shared_lib()
        install_.run(self)


if __name__ == "__main__":
    cmdclass = {}
    try:
        from wheel.bdist_wheel import bdist_wheel as bdist_wheel_
        class bdist_wheel(bdist_wheel_):
            def run(self):
                prepare_shared_lib()
                bdist_wheel_.run(self)

            def finalize_options(self):
                bdist_wheel_.finalize_options(self)
                self.universal = True
                self.plat_name_supplied = True
                self.plat_name = get_platform()
        cmdclass = {'bdist_wheel': bdist_wheel, 'install': install}
    except ImportError:
        cmdclass = {'install': install}

    # cibuildwheel for windows wasn't properly cleaning up after each build
    # it is necessary to fully remove the build directory once at the beginning
    # prior to building.
    if platform.system() == "Windows":
        if os.path.exists("build"):
            print("Removing old build directory")
            shutil.rmtree("build")
        else:
            print("No previous build directory detected")

    long_description = ""
    try:
        long_description = open("README.rst", encoding="utf8").read()
    except:
        long_description = open("README.rst").read()
    setup(
        name=NAME,
        description=find_meta("description"),
        long_description=long_description,
        license=find_meta("license"),
        url=URI,
        version=find_meta("version"),
        author=find_meta("author"),
        author_email=find_meta("email"),
        maintainer=find_meta("author"),
        maintainer_email=find_meta("email"),
        keywords=KEYWORDS,
        packages=PACKAGES,
        package_data={"parasail": [get_libname()] + get_includes()},
        cmdclass=cmdclass,
        zip_safe=False,
        classifiers=CLASSIFIERS,
        install_requires=INSTALL_REQUIRES,
    )

