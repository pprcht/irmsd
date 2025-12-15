import bz2
import gzip
import lzma
import os
import pickle

import pytest

pytest.importorskip("ase")

from irmsd.utils.io import dump_results_to_pickle, read_structures

CAFFEINE_SINGLE_SDF = """
  Avogadro

 24 25  0  0  0  0  0  0  0  0999 V2000
    1.1426    0.9750   -0.1011 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.3032   -0.4158   -0.1836 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1503   -1.1865   -0.1294 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0770    1.6825    0.0339 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1838    0.7902    0.0760 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0689   -0.5818   -0.0002 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5402    1.0307    0.2021 N   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1379   -0.1838    0.1918 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2791   -1.1896    0.0708 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.2860   -2.6331   -0.2127 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4056   -0.9268   -0.2977 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1247    2.9069    0.1020 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1825    2.3266    0.3154 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3788    1.7466   -0.1583 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2060   -0.2943    0.2751 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7647   -2.9107   -1.1477 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7056   -3.0647   -0.1590 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.8977   -3.0000    0.6068 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9782    2.9276   -0.5655 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8075    2.8576    1.1847 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2527    2.1730    0.4131 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.8966    1.5565   -1.0942 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.0359    1.4620    0.6585 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.1148    2.7928   -0.0810 H   0  0  0  0  0  0  0  0  0  0  0  0
 10 16  1  0  0  0  0
 14 22  1  0  0  0  0
 13 19  1  0  0  0  0
  2 11  2  0  0  0  0
 10 17  1  0  0  0  0
  3 10  1  0  0  0  0
 10 18  1  0  0  0  0
  2  3  1  0  0  0  0
  1  2  1  0  0  0  0
  1 14  1  0  0  0  0
 14 24  1  0  0  0  0
 14 23  1  0  0  0  0
  3  6  1  0  0  0  0
  1  4  1  0  0  0  0
  6  9  1  0  0  0  0
  5  6  2  0  0  0  0
  4  5  1  0  0  0  0
  4 12  2  0  0  0  0
  8  9  2  0  0  0  0
  5  7  1  0  0  0  0
  7  8  1  0  0  0  0
  8 15  1  0  0  0  0
  7 13  1  0  0  0  0
 13 21  1  0  0  0  0
 13 20  1  0  0  0  0
M  END
> <modelView>


> <projection>


$$$$
"""

PEPTIDE_CONFORMER_SDF = """structure 1 of 4
 crest11282510233D Energy =       -3.29043291

 20 19  0  0  0  0  0  0  0  0999 V2000
   -1.6320   -0.4951   -0.4391 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0605    0.7861    0.1356 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1685    1.2698   -0.2693 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.1632    0.5066   -1.0143 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0698   -0.2507   -0.0882 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7122   -1.3113    0.3694 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1712   -1.6261    0.3533 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7104    1.3592    0.9871 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2566    0.3186    0.2031 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4965   -1.5112    1.3023 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1501   -1.6022    0.3892 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1435   -0.4339   -0.4795 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2613   -0.6285   -1.4667 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.5690    2.0387    0.2655 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.7599    1.1927   -1.6228 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6795   -0.1841   -1.7062 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.7216   -0.2775    0.8201 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4920    0.3764   -1.1206 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5684   -1.3722   -0.8364 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5306   -0.2456    0.5218 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  7  1  0  0  0  0
  1 12  1  0  0  0  0
  1 13  1  0  0  0  0
  2  3  1  0  0  0  0
  2  8  2  0  0  0  0
  3  4  1  0  0  0  0
  3 14  1  0  0  0  0
  4  5  1  0  0  0  0
  4 15  1  0  0  0  0
  4 16  1  0  0  0  0
  5  6  2  0  0  0  0
  5  9  1  0  0  0  0
  7 10  1  0  0  0  0
  7 11  1  0  0  0  0
  9 17  1  0  0  0  0
 12 18  1  0  0  0  0
 12 19  1  0  0  0  0
 12 20  1  0  0  0  0
M  END
$$$$
structure 2 of 4
 crest11282510233D Energy =       -3.29029914

 20 19  0  0  0  0  0  0  0  0999 V2000
   -1.9260   -0.5887   -0.4370 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0374    0.5997   -0.1446 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0848    0.7532   -0.9410 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.3162    1.1783   -0.3121 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1261   -0.0283    0.0798 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6054   -1.1203    0.1279 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3471   -1.7410    0.2513 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2674    1.3013    0.8193 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4258    0.1822    0.3706 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5528   -1.6656    1.2374 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3287   -1.6557    0.1821 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3579   -0.3296   -0.0286 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8839   -0.7972   -1.5182 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.2734    0.0410   -1.6479 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.0619    1.7643    0.5763 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.8850    1.8215   -0.9895 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.8138   -0.6804    0.6108 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7862    0.5013   -0.5898 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9794   -1.2120   -0.1813 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3993   -0.0636    1.0274 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  7  1  0  0  0  0
  1 12  1  0  0  0  0
  1 13  1  0  0  0  0
  2  3  1  0  0  0  0
  2  8  2  0  0  0  0
  3  4  1  0  0  0  0
  3 14  1  0  0  0  0
  4  5  1  0  0  0  0
  4 15  1  0  0  0  0
  4 16  1  0  0  0  0
  5  6  2  0  0  0  0
  5  9  1  0  0  0  0
  7 10  1  0  0  0  0
  7 11  1  0  0  0  0
  9 17  1  0  0  0  0
 12 18  1  0  0  0  0
 12 19  1  0  0  0  0
 12 20  1  0  0  0  0
M  END
$$$$
structure 3 of 4
 crest11282510233D Energy =       -3.28845436

 20 19  0  0  0  0  0  0  0  0999 V2000
   -2.0478   -0.4510   -0.3834 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0699    0.6977   -0.2281 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0612    0.6730   -1.0242 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.3175    1.1076   -0.4522 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0251   -0.0576    0.1871 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4651   -1.1244    0.3054 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2904   -1.6566   -0.7382 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2623    1.5571    0.6089 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2865    0.1618    0.6101 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7316   -2.4981   -0.3940 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3499   -1.5805   -0.3289 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8454   -0.6335    0.8895 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7234   -0.2192   -1.2221 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.1906   -0.1144   -1.6693 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.1012    1.8704    0.3020 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.9460    1.5650   -1.2215 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.6098   -0.6703    1.0044 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1774   -0.8843    1.7130 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3588    0.2901    1.1597 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5843   -1.4292    0.7952 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  7  1  0  0  0  0
  1 12  1  0  0  0  0
  1 13  1  0  0  0  0
  2  3  1  0  0  0  0
  2  8  2  0  0  0  0
  3  4  1  0  0  0  0
  3 14  1  0  0  0  0
  4  5  1  0  0  0  0
  4 15  1  0  0  0  0
  4 16  1  0  0  0  0
  5  6  2  0  0  0  0
  5  9  1  0  0  0  0
  7 10  1  0  0  0  0
  7 11  1  0  0  0  0
  9 17  1  0  0  0  0
 12 18  1  0  0  0  0
 12 19  1  0  0  0  0
 12 20  1  0  0  0  0
M  END
$$$$
structure 4 of 4
 crest11282510233D Energy =       -3.28768690

 20 19  0  0  0  0  0  0  0  0999 V2000
   -1.2329   -0.6937    0.1791 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9849    0.7929    0.3066 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2562    1.2903    0.6382 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.4360    0.4791    0.9175 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9105   -0.2495   -0.3052 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3658   -1.2776   -0.6344 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2224   -1.0960    1.1650 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9326    1.5251    0.1089 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9296    0.3037   -0.9931 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0645   -0.5652    0.9904 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9113   -0.8416    2.0940 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6734   -1.0307   -1.2278 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3065   -1.2470    0.3935 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.3953    2.3002    0.6303 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.2019   -0.2350    1.7106 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.2449    1.1140    1.2911 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.1157   -0.2755   -1.7565 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5537   -0.4457   -1.4938 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8867   -0.8011   -1.9482 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9257   -2.0874   -1.3171 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  7  1  0  0  0  0
  1 12  1  0  0  0  0
  1 13  1  0  0  0  0
  2  3  1  0  0  0  0
  2  8  2  0  0  0  0
  3  4  1  0  0  0  0
  3 14  1  0  0  0  0
  4  5  1  0  0  0  0
  4 15  1  0  0  0  0
  4 16  1  0  0  0  0
  5  6  2  0  0  0  0
  5  9  1  0  0  0  0
  7 10  1  0  0  0  0
  7 11  1  0  0  0  0
  9 17  1  0  0  0  0
 12 18  1  0  0  0  0
 12 19  1  0  0  0  0
 12 20  1  0  0  0  0
M  END
$$$$
"""


@pytest.fixture(scope="module")
def tmp_dir_fixture(tmp_path_factory):
    """Fixture that provides a temporary directory for tests."""
    tmp_path = tmp_path_factory.mktemp("io_test_files")
    return tmp_path


@pytest.fixture(scope="module")
def caffeine_sdf_file_fixture(tmp_dir_fixture):
    """Fixture that provides a temporary SDF file for caffeine."""
    sdf_path = tmp_dir_fixture / "caffeine.sdf"

    with open(sdf_path, "w") as f:
        f.write(CAFFEINE_SINGLE_SDF)
    return str(sdf_path)


@pytest.fixture(scope="module")
def peptide_conformer_sdf_file_fixture(tmp_dir_fixture):
    """Fixture that provides a temporary SDF file for peptide conformers."""
    sdf_path = tmp_dir_fixture / "peptide_conformers.sdf"
    with open(sdf_path, "w") as f:
        f.write(PEPTIDE_CONFORMER_SDF)
    return str(sdf_path)


# xyz file types are already testet via other tests such as cli tests
@pytest.mark.parametrize(
    "file_fixture",
    [
        "caffeine_sdf_file_fixture",
        "peptide_conformer_sdf_file_fixture",
    ],
)
def test_read_structures_ase(file_fixture, request):
    """Test the read_structures function with ASE-supported files."""
    file = request.getfixturevalue(file_fixture)
    structures = read_structures([file])
    assert len(structures) > 0
    assert all(structure.get_atomic_numbers() is not None for structure in structures)


@pytest.fixture
def molecules():
    # simple pickleable stand-ins
    return [{"id": 1}, {"id": 2}]


@pytest.fixture
def results():
    return {"energy": -123.4, "status": "ok"}


def load_pickle(path, compress):
    """Helper to load a pickle with optional compression."""
    if compress is None:
        open_fn = open
    elif compress == "gz":
        open_fn = gzip.open
    elif compress == "bz2":
        open_fn = bz2.open
    elif compress == "xz":
        open_fn = lzma.open
    else:
        raise AssertionError("unexpected compress")

    with open_fn(path, "rb") as f:
        return pickle.load(f)


def test_dump_without_results_no_compression(tmp_path, molecules):
    outfile = tmp_path / "test.pkl"

    final_path = dump_results_to_pickle(
        molecules=molecules,
        outfile=str(outfile),
        results=None,
        compress=None,
    )

    assert final_path == str(outfile)
    assert os.path.exists(final_path)

    payload = load_pickle(final_path, compress=None)
    assert payload == molecules


def test_dump_with_results_no_compression(tmp_path, molecules, results):
    outfile = tmp_path / "test.pkl"

    final_path = dump_results_to_pickle(
        molecules=molecules,
        outfile=str(outfile),
        results=results,
        compress=None,
    )

    payload = load_pickle(final_path, compress=None)

    assert isinstance(payload, dict)
    assert payload["molecules"] == molecules
    for k, v in results.items():
        assert payload[k] == v


@pytest.mark.parametrize("compress", ["gz", "bz2", "xz"])
def test_dump_with_compression(tmp_path, molecules, results, compress):
    outfile = tmp_path / "test.pkl"

    final_path = dump_results_to_pickle(
        molecules=molecules,
        outfile=str(outfile),
        results=results,
        compress=compress,
    )

    assert final_path.endswith(f".pkl.{compress}")
    assert os.path.exists(final_path)

    payload = load_pickle(final_path, compress=compress)
    assert payload["molecules"] == molecules
    assert payload["energy"] == results["energy"]


def test_outfile_extension_is_normalized(tmp_path, molecules):
    # wrong extension should be replaced by .pkl
    outfile = tmp_path / "data.txt"

    final_path = dump_results_to_pickle(
        molecules=molecules,
        outfile=str(outfile),
        compress=None,
    )

    assert final_path.endswith(".pkl")
    assert os.path.exists(final_path)


def test_double_extension_handling(tmp_path, molecules):
    # foo.pkl.gz with compress=None should normalize to foo.pkl
    outfile = tmp_path / "foo.pkl.gz"

    final_path = dump_results_to_pickle(
        molecules=molecules,
        outfile=str(outfile),
        compress=None,
    )

    assert final_path.endswith("foo.pkl")
    assert os.path.exists(final_path)


def test_results_with_molecules_key_raises(tmp_path, molecules):
    outfile = tmp_path / "test.pkl"
    bad_results = {"molecules": "oops"}

    with pytest.raises(ValueError, match="contains a key 'molecules'"):
        dump_results_to_pickle(
            molecules=molecules,
            outfile=str(outfile),
            results=bad_results,
        )


def test_invalid_compression_raises(molecules, tmp_path):
    outfile = tmp_path / "test.pkl"

    with pytest.raises(ValueError, match="Invalid compress"):
        dump_results_to_pickle(
            molecules=molecules,
            outfile=str(outfile),
            compress="zip",
        )


def test_existing_pkl_extension_with_compression(tmp_path, molecules):
    outfile = tmp_path / "data.pkl"

    final_path = dump_results_to_pickle(
        molecules=molecules,
        outfile=str(outfile),
        compress="gz",
    )

    assert final_path.endswith(".pkl.gz")
    payload = load_pickle(final_path, compress="gz")
    assert payload == molecules
