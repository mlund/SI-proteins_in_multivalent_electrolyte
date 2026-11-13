from analysis import run


def test_runs_analysis_notebooks_but_never_the_mc_simulation():
    names = [p.name for p in run.analysis_notebooks()]
    assert names
    assert all(p.exists() for p in run.analysis_notebooks())
    assert not any("simulate" in n for n in names)
