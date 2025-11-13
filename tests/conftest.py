import numpy as np
import pytest

CAFFEINE_XTB = """24

C          1.07317        0.04885       -0.07573
N          2.51365        0.01256       -0.07580
C          3.35199        1.09592       -0.07533
N          4.61898        0.73028       -0.07549
C          4.57907       -0.63144       -0.07531
C          3.30131       -1.10256       -0.07524
C          2.98068       -2.48687       -0.07377
O          1.82530       -2.90038       -0.07577
N          4.11440       -3.30433       -0.06936
C          5.45174       -2.85618       -0.07235
O          6.38934       -3.65965       -0.07232
N          5.66240       -1.47682       -0.07487
C          7.00947       -0.93648       -0.07524
C          3.92063       -4.74093       -0.06158
H          0.73398        1.08786       -0.07503
H          0.71239       -0.45698        0.82335
H          0.71240       -0.45580       -0.97549
H          2.99301        2.11762       -0.07478
H          7.76531       -1.72634       -0.07591
H          7.14864       -0.32182        0.81969
H          7.14802       -0.32076       -0.96953
H          2.86501       -5.02316       -0.05833
H          4.40233       -5.15920        0.82837
H          4.40017       -5.16929       -0.94780
"""

CAFFEINE_XTB_ROTATED = """24

C      0.000000    0.000000    0.000000
N      0.000000   -0.000000    1.440937
C      0.027298    1.103793    2.251727
N      0.018994    0.770281    3.527523
C     -0.015495   -0.591576    3.521921
C     -0.028175   -1.094567    2.256431
C     -0.064687   -2.486034    1.970767
O     -0.073776   -2.928423    0.826168
N     -0.089014   -3.274320    3.124715
C     -0.073958   -2.792858    4.450344
O     -0.093669   -3.572210    5.407882
N     -0.036585   -1.409132    4.626198
C     -0.021822   -0.835227    5.959233
C     -0.133076   -4.714691    2.967187
H      0.025262    1.029828   -0.365250
H     -0.911742   -0.491950   -0.347970
H      0.886557   -0.536074   -0.347902
H      0.052263    2.115820    1.867129
H     -0.040596   -1.605572    6.734725
H     -0.900910   -0.194914    6.082835
H      0.887769   -0.238932    6.082275
H     -0.144049   -5.023239    1.919009
H     -1.032993   -5.098154    3.459225
H      0.742359   -5.153025    3.457406
"""

CAFFEINE_OBABEL = """24

C      0.8386861685    2.3213832518    0.0010145738
N      1.0653014660    0.8729242511   -0.0000214594
C      2.3184263459    0.3223786159   -0.0008480483
N      1.6274263583   -0.5517078054   -0.0022388503
C      0.7215992891   -1.2013751523   -0.0019368855
C      0.0362179478   -0.0158207422    0.0003381297
C     -1.4349785264    0.0369662984    0.0033774333
O     -2.0569094114    1.0837423228    0.0062761005
N     -2.0392433349   -1.1889082038    0.0038185877
C     -1.3887949619   -2.3933182627    0.0004337555
O     -2.0022408409   -3.4463203653    0.0006547338
N     -0.0244704068   -2.3686449774   -0.0021002133
C      0.7646931543   -3.6133446731   -0.0036790374
C     -3.5211821920   -1.2424968568    0.0094433313
H      1.7927742350    2.8549060126    0.0032623100
H      0.2677844604    2.5994076569    0.8913082227
H      0.2711447399    2.6009557963   -0.8909346942
H      3.2197333040    0.9168049818    0.0008472191
H      0.1266508985   -4.5003129988   -0.0062220724
H      1.3971083521   -3.6365855905    0.8877087428
H      1.3994336185   -3.6333386595   -0.8933820058
H     -3.9752819067   -0.2490883471    0.0122296564
H     -3.8576064983   -1.7785508485    0.9010397951
H     -3.8648098649   -1.7768842904   -0.8804701717
"""

# FIXME: Check with crest whether this is the correct structure
CAFFEINE_OBABEL_ALIGNED = """24

C      1.154606     0.109937    -0.076297
N      2.620590     0.093227    -0.075995
C      3.371311     1.237709    -0.076716
N      4.118713     0.411331    -0.077013
C      4.609179    -0.589685    -0.075762
C      3.326384    -1.068984    -0.074410
C      3.030366    -2.511057    -0.070913
O      1.894949    -2.950797    -0.068823
N      4.139067    -3.310234    -0.069063
C      5.434664    -2.868512    -0.071495
O      6.371364    -3.648077    -0.070031
N      5.636572    -1.518986    -0.074525
C      6.994902    -0.947149    -0.075159
C      3.946168    -4.780540    -0.062873
H      0.786678     1.139288    -0.074902
H      0.784954    -0.406510     0.813921
H      0.785602    -0.403837    -0.968323
H      2.934571     2.225108    -0.075914
H      7.763790    -1.723438    -0.076614
H      7.121882    -0.326894     0.816032
H      7.120681    -0.324959    -0.965062
H      2.891210    -5.063621    -0.060902
H      4.418204    -5.200748     0.829363
H      4.416983    -5.208472    -0.952145
"""


@pytest.fixture(scope="session")
def caffeine_molecule_xyz():
    """Fixture that returns a caffeine molecule."""
    return CAFFEINE_XTB


@pytest.fixture(scope="session")
def caffeine_molecule_xyz_rotated():
    """Fixture that returns a caffeine molecule."""
    return CAFFEINE_XTB_ROTATED


@pytest.fixture(scope="session")
def caffeine_molecule_xyz_obabel():
    """Fixture that returns a caffeine molecule."""
    return CAFFEINE_OBABEL


@pytest.fixture(
    scope="session",
    params=[
        # Conformer1, Conformer2, mask, rmsd, aligned_conformer, Umat
        (
            CAFFEINE_XTB,
            CAFFEINE_XTB_ROTATED,
            0.00000042,
            CAFFEINE_XTB,
            # FIXME: check with crest, this is actually taken from this code itself but the rmsd is the same as crest
            np.asarray(
                [
                    [0.000586, 0.025178, 0.999683],
                    [0.025178, 0.999366, -0.025185],
                    [-0.999683, 0.025185, -0.000049],
                ],
                dtype=np.float64,
            ),
        ),
        (
            CAFFEINE_XTB,
            CAFFEINE_OBABEL,
            0.13939435,
            CAFFEINE_OBABEL_ALIGNED,
            # FIXME: check with crest, this is actually taken from this code itself but the rmsd is the same as crest
            np.asarray(
                [
                    [0.165823, -0.986155, -0.000907],
                    [0.986155, 0.165823, 0.000503],
                    [-0.000346, -0.000978, 0.999999],
                ],
                dtype=np.float64,
            ),
        ),
    ],
)
def caffeine_rmsd_test_data(request):
    """Fixture that returns atom numbers and positions for caffeine molecule in
    different conformations."""
    return request.param


@pytest.fixture(
    scope="session",
    params=[
        # Conformer1, rot constants in MHz, avmom in a.u., evec
        (
            CAFFEINE_XTB,
            # FIXME: check with crest, this is actually taken from this code itself but the rmsd is the same as crest
            np.asarray(
                [
                    1068.073127,
                    710.711798,
                    430.26121,
                ],
                dtype=np.float64,
            ),
            1.30564546e-44,
            np.asarray(
                [
                    [0.671158, -0.741313, 0.001206],
                    [0.741314, 0.671157, -0.001149],
                    [0.000042, 0.001665, 0.999999],
                ],
                dtype=np.float64,
            ),
        ),
    ],
)
def caffeine_axis_test_data(request):
    return request.param


@pytest.fixture(
    scope="session",
    params=[
        # Conformer1, cn
        (
            CAFFEINE_OBABEL,
            # FIXME: check with crest, this is actually taken from this code itself but the rmsd is the same as crest
            np.asarray(
                # fmt: off
            [3.689705,3.866219,2.934731,3.626334,3.144874,3.791101,2.768159,0.858364,2.743387,
             2.716063,0.860503,2.747386,3.698299,3.698564,0.924137,0.924083,0.924084,0.925654,
             0.924202,0.92413,0.924141,0.924245,0.924105,0.924098],
                # fmt: on
                dtype=np.float64,
            ),
        ),
    ],
)
def caffeine_cn_test_data(request):
    return request.param


@pytest.fixture(
    scope="session",
    params=[
        # Conformer1, heavy, canonical_atom_id
        (
            CAFFEINE_XTB,
            False,
            # FIXME: check with crest, this is actually taken from this code itself but the rmsd is the same as crest
            np.asarray(
                # fmt: off
                [1,8,5,7,13,14,12,6,10,9,2,11,4,3,15,15,15,18,17,17,17,16,16,16],
                # fmt: on
                dtype=int,
            ),
        ),
        (
            CAFFEINE_XTB,
            True,
            # FIXME: check with crest, this is actually taken from this code itself but the rmsd is the same as crest
            np.asarray(
                # fmt: off
                [1,8,5,7,13,14,12,6,10,9,2,11,4,3,0,0,0,0,0,0,0,0,0,0],
                # fmt: on
                dtype=int,
            ),
        ),
    ],
)
def caffeine_canonical_test_data(request):
    return request.param
