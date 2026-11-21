"""Run the two semantic controls on the same captured pack, offline."""

from pathlib import Path
import subprocess
import sys
import tempfile


def run(args, expected_code, markers):
    cmd = [sys.executable, "-B", str(Path(__file__).with_name("check.py")), *args]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, check=False)
    if result.returncode != expected_code or any(s not in result.stdout for s in markers):
        raise RuntimeError(f"unexpected checker outcome ({result.returncode})\n{result.stdout}")
    return result.stdout


def main():
    clean = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).with_name("pack.json.gz")
    transitions = ["integrity PASS", "all local transitions PASS"]
    with tempfile.TemporaryDirectory() as directory:
        changed = Path(directory) / "double-exit-gain.json.gz"
        run([clean, "--endpoint-only"], 0, transitions + ["completed-pack Close PASS"])
        print("PASS clean tape: rehash, every recorded transition, padded complex mean")
        changed_id = run([clean, "--negative", changed], 0, ["generated semantic negative"])
        print(changed_id.strip())
        output = run([changed, "--endpoint-only"], 1, transitions + ["REJECT padded ternary mean"])
        print(next(line for line in output.splitlines() if line.startswith("REJECT")))
        print("PASS semantic negative: fresh hash and consistent rounded states; wrong mean")


if __name__ == "__main__":
    main()
