import subprocess
import unittest


class TestCliSmoke(unittest.TestCase):
    def run_cmd(self, cmd):
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if proc.returncode != 0:
            self.fail(
                "Command failed: {cmd}\nexit={code}\nstdout:\n{out}\nstderr:\n{err}".format(
                    cmd=" ".join(cmd),
                    code=proc.returncode,
                    out=proc.stdout,
                    err=proc.stderr,
                )
            )

    def test_help_commands(self):
        self.run_cmd(["python3", "compute_rates_misfit.py", "--help"])
        self.run_cmd(["python3", "compute_rates_single.py", "--help"])
        self.run_cmd(["python3", "create_trench_motion_table.py", "--help"])

    def test_make_smoke(self):
        self.run_cmd(["make", "smoke"])

    def test_make_single_smoke(self):
        self.run_cmd(["make", "single-smoke"])


if __name__ == "__main__":
    unittest.main()
