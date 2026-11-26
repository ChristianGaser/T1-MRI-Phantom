"""Download of the model files t1sim needs.

The front end segments a real scan with CAT-Surface's unified segmentation,
whose tissue priors are a series of five shooting templates with six classes
(GM, WM, CSF, bone, soft tissue, air).  The templates belong to T1Prep's
source tree but are left out of its PyPI package, so t1sim fetches them
itself, from the commit of the T1Prep repository that added them, and checks
every file against its SHA-256.

Where they go (:func:`model_dir`):

- ``T1SIM_MODEL_DIR`` when set;
- the folder ``models`` of this repository when t1sim runs from a source
  checkout;
- the user cache otherwise (an installed package has no repository, and
  files written into a package folder are left behind by pip).

``python -m t1sim.models`` downloads what is missing; ``--verify`` hashes
the files that are there as well.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import sys
from pathlib import Path
from urllib.parse import urlsplit
from urllib.request import urlopen

#: The commit of the T1Prep repository the templates are taken from.
MODEL_SOURCE_COMMIT = "2635ddcfd2d47a5f4623657a5c7ff960c80a3b8e"

_DEFAULT_MODEL_URL = (
    "https://raw.githubusercontent.com/ChristianGaser/T1Prep/"
    f"{MODEL_SOURCE_COMMIT}/src/t1prep/data/templates_MNI152NLin2009cAsym"
)

#: The files and their SHA-256: the six-class shooting templates, coarse to fine.
MODEL_FILES = {
    "Template_0_GS_6class.nii.gz": "e2898ee33f6068b2e09b0eb2b55991e2fcff33e5fdf6c5aaa29e1f50cc79c999",
    "Template_1_GS_6class.nii.gz": "fafdae2cae27688748a0b0ae216c54b54e1ca0cd872e04865c72dd8ec0627cb6",
    "Template_2_GS_6class.nii.gz": "07bb84cc11e41e2aedbc5074fed7aa74c59f5a9e01e582767231565a6b1dd9df",
    "Template_3_GS_6class.nii.gz": "454c698c985ab85af33d3ea19664e84e6fecea0711d924ef819cfec9838d4f86",
    "Template_4_GS_6class.nii.gz": "bccf3f8109d98aa4a43c06e82283b70c36389d05f4448b208986b299b8457586",
}


def model_url() -> str:
    """Folder the files are downloaded from; ``T1SIM_MODEL_URL`` overrides it."""
    return os.environ.get("T1SIM_MODEL_URL", _DEFAULT_MODEL_URL).rstrip("/")


def _repository_root():
    """The root of the source checkout this file lies in, or None for an installed package."""
    root = Path(__file__).resolve().parents[2]
    if (root / "pyproject.toml").is_file() and (root / "src" / "t1sim").is_dir():
        return root
    return None


def _user_cache_dir() -> Path:
    if sys.platform == "win32":
        base = os.environ.get("LOCALAPPDATA")
        return (Path(base) if base else Path.home() / "AppData" / "Local") / "t1sim" / "Cache"
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Caches" / "t1sim"
    xdg = os.environ.get("XDG_CACHE_HOME")
    return Path(xdg) / "t1sim" if xdg else Path.home() / ".cache" / "t1sim"


def model_dir() -> Path:
    """The folder that holds the model files (see the module notes)."""
    override = os.environ.get("T1SIM_MODEL_DIR")
    if override:
        return Path(override)
    root = _repository_root()
    if root is not None:
        return root / "models"
    return _user_cache_dir() / "models"


def sha256(path) -> str:
    """SHA-256 of a file as a hex string."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def missing_models(directory=None, verify: bool = False) -> list:
    """Names of the model files that are not in *directory*.

    Args:
        directory: Folder to look in; default: :func:`model_dir`.
        verify: Also count a file whose SHA-256 differs as missing.  Hashes
            all files (44 MB), so it is off for the check before every run.
    """
    directory = Path(directory) if directory is not None else model_dir()
    missing = []
    for name, digest in MODEL_FILES.items():
        path = directory / name
        if not path.is_file() or (verify and sha256(path) != digest):
            missing.append(name)
    return missing


def _fetch(url: str, dest: Path) -> None:
    """Write the content of *url* to *dest*.

    http(s) goes through requests when it is installed: it verifies against
    certifi's CA bundle, which python.org's macOS Python needs (its OpenSSL
    has none until "Install Certificates.command" has run, and urllib then
    fails with CERTIFICATE_VERIFY_FAILED).
    """
    if urlsplit(url).scheme.lower() in ("http", "https"):
        try:
            import requests
        except ImportError:
            requests = None
        if requests is not None:
            with requests.get(url, stream=True, timeout=60) as response:
                response.raise_for_status()
                with open(dest, "wb") as out:
                    for chunk in response.iter_content(chunk_size=1 << 20):
                        out.write(chunk)
            return
    with urlopen(url, timeout=60) as response, open(dest, "wb") as out:
        shutil.copyfileobj(response, out)


def download_models(directory=None, force: bool = False, verbose: bool = True) -> Path:
    """Download the model files and check them against their SHA-256.

    A file is written under a temporary name and renamed once its hash is
    right, so an interrupted or corrupted download never passes for a model.

    Args:
        directory: Folder to download into; default: :func:`model_dir`.
        force: Download the files that are there as well.
        verbose: Print what is downloaded.

    Returns:
        The folder with the files.

    Raises:
        RuntimeError: A download failed or a file has the wrong hash.
    """
    directory = Path(directory) if directory is not None else model_dir()
    directory.mkdir(parents=True, exist_ok=True)
    base = model_url()
    names = list(MODEL_FILES) if force else missing_models(directory)
    if names and verbose:
        print(f"Downloading {len(names)} model file(s) to {directory}\n  from {base}")
    for name in names:
        url = f"{base}/{name}"
        # Named after the process, so two runs that start together do not
        # write into one file; opened as any file is, so it gets the usual
        # permissions (mkstemp would leave the model readable by its owner only)
        tmp = directory / f"{name}.{os.getpid()}.part"
        try:
            try:
                _fetch(url, tmp)
            except Exception as exc:
                raise RuntimeError(
                    f"Failed to download {url}: {exc}\n"
                    "Check the internet connection, or set T1SIM_MODEL_URL to "
                    "another folder that holds the files (file:// works)."
                ) from exc
            digest = sha256(tmp)
            if digest != MODEL_FILES[name]:
                raise RuntimeError(
                    f"{name} from {url} has SHA-256 {digest}, "
                    f"expected {MODEL_FILES[name]}: the file was not kept."
                )
            os.replace(tmp, directory / name)
        finally:
            if tmp.exists():
                tmp.unlink()
        if verbose:
            size = (directory / name).stat().st_size / 1e6
            print(f"  {name} ({size:.1f} MB)")
    return directory


def ensure_models(directory=None, verbose: bool = True) -> Path:
    """The folder with the model files, downloading those that are missing.

    Args:
        directory: Folder to use; default: :func:`model_dir`.
        verbose: Print what is downloaded.
    """
    directory = Path(directory) if directory is not None else model_dir()
    if missing_models(directory):
        download_models(directory, verbose=verbose)
    return directory


def main(argv=None) -> int:
    """Command line of ``python -m t1sim.models``."""
    parser = argparse.ArgumentParser(
        prog="python -m t1sim.models",
        description="Download the model files of t1sim (the six-class shooting "
                    "templates of the unified segmentation, 44 MB).",
    )
    parser.add_argument("--dir", default=None,
                        help="folder to download into (default: see --where)")
    parser.add_argument("--force", action="store_true",
                        help="download the files that are there as well")
    parser.add_argument("--verify", action="store_true",
                        help="hash the files that are there and download those that differ")
    parser.add_argument("--where", action="store_true",
                        help="print the model folder and exit")
    args = parser.parse_args(argv)
    directory = Path(args.dir) if args.dir else model_dir()
    if args.where:
        print(directory)
        return 0
    if args.verify:
        for name in missing_models(directory, verify=True):
            path = directory / name
            if path.is_file():
                print(f"{name}: wrong SHA-256, downloading it again")
                path.unlink()
    try:
        download_models(directory, force=args.force)
    except RuntimeError as exc:
        print(exc, file=sys.stderr)
        return 1
    print(f"All {len(MODEL_FILES)} model files are in {directory}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
