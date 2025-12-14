import subprocess
import filecmp
import os
import pytest
from pyinstrument import Profiler


@pytest.fixture
def pyinstrument_profile_html(request):
    profiler = Profiler()
    profiler.start()
    yield profiler
    profiler.stop()
    html_report = profiler.output_html()
    report_path = os.path.join(os.path.dirname(__file__), "pyinstrument_report.html")
    with open(report_path, "w") as f:
        f.write(html_report)
    print(f"\nPyInstrument HTML report written to: {report_path}")


def test_output_main_tsv_regression(tmp_path, pyinstrument_profile_html):
    # Run the main script
    subprocess.run(["python", "main.py"], check=True, cwd=os.path.dirname(__file__))

    # Compare against the snapshot
    expected = os.path.join(os.path.dirname(__file__), "snapshot", "output.tsv")
    actual = os.path.join(os.path.dirname(__file__), "output_main.tsv")

    assert os.path.exists(actual), "output_main.tsv was not created"
    assert filecmp.cmp(expected, actual, shallow=False), (
        "output_main.tsv does not match snapshot/output.tsv"
    )
