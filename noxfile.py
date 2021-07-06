import nox


def setup_test_env(session):
    session.install("--upgrade", "pip", "setuptools")
    session.install("-r", "requirements.txt", "-r", "requirements-dev.txt")
    session.run("pip", "install", "--no-deps", "--editable", ".")


@nox.session(python=["3.6", "3.7", "3.8", "3.9"])
def test(session):
    setup_test_env(session)
    session.run("pytest", "--doctest-modules", "test/")


@nox.session(python="3.9")
def black(session):
    session.install("black")
    session.run("black", "--check", "--diff", "scripts/", "test/")


@nox.session
def lint(session):
    setup_test_env(session)
    session.run("flake8", "scripts", "test")
