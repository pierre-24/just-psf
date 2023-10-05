import pathlib
import shutil


def path_from_tests_files(path: pathlib.Path) -> pathlib.Path:
    p = pathlib.Path(__file__).parent / path

    assert p.exists()
    return p


def copy_to_temporary_directory(source: pathlib.Path, path_tempdir: pathlib.Path) -> pathlib.Path:
    """Copy the content of a file to the temporary directory
    """

    full_path_source = path_from_tests_files(source)

    path_dest = path_tempdir / full_path_source.name
    assert not path_dest.exists()

    shutil.copy(full_path_source, path_dest)
    return path_dest
