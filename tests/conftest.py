# import fixtures that setup data for the tests
# solution taken from https://gist.github.com/peterhurford/09f7dcda0ab04b95c026c60fa49c2a68
from glob import glob


def _as_module(fixture_path: str) -> str:
    return fixture_path.replace("/", ".").replace("\\", ".").replace(".py", "")


pytest_plugins = [_as_module(fixture) for fixture in glob("tests/fixtures/[!_]*.py")]
