"""Observe the pinned diagnostic through the return of its left forward transform."""

import argparse
from collections import Counter
from contextlib import ExitStack
from fractions import Fraction
import gzip
import hashlib
from io import BytesIO
import json
from pathlib import Path
import sys
from unittest.mock import patch


COMMIT = "9585c1fa2ed5f2316ac09e04b357cb720badfcb8"
PACK = "3x5-P72/left-forward/pack-1/physical-axis-1-side-3"


class PackComplete(Exception):
    """Leave the published call after the entire observed transform has returned."""


def require(condition, message):
    if not condition:
        raise ValueError(message)


def snapshot(tape, width, total):
    # Bypass the meter for observation only. Always restore the cursor.
    stream = getattr(tape, "stream", tape)
    saved = stream.tell()
    try:
        stream.seek(0)
        data = stream.read()
    finally:
        stream.seek(saved)
    require(len(data) == 2 * width * total, "unexpected physical tape length")
    values = [
        int.from_bytes(data[i:i + width], "little", signed=True)
        for i in range(0, len(data), width)
    ]
    return list(zip(values[::2], values[1::2], strict=True))


def observe(source):
    sub = source / "research" / "triangular-native-supplier"
    sys.path.insert(0, str(sub.resolve()))
    import multiplier
    import product_path
    import measure_product

    package = product_path.fourier()
    report = dict(commit=COMMIT, pack_id=PACK, diagnostic="measure_product.measure(Product,72)",
                  axis=1, pack_index=1, operand="left", gates=[], states=[], tables=[])
    origin = {}
    tables = {}
    kernel_calls = Counter()
    axis_calls = []
    wrap_product = product_path.Product.multiply_tapes
    wrap_axis = package.fourier_axis
    spectral = sys.modules[wrap_axis.__module__]
    wrap_polar = spectral.polar_root
    wrap_word = package.Grid.word
    wrap_diag = multiplier.tensor_diagonal
    wrap_h = product_path.Product.h
    wrap_transform = multiplier.Product.transform
    wrap_crt = product_path.crt_tape
    wrap_hadamard = product_path.hadamard_tape
    machine = None
    selected = False
    finished = False

    def polar_root(compiler, radii, turns, bits, phase_grid=None):
        root, phases = wrap_polar(compiler, radii, turns, bits, phase_grid)
        if len(turns) == 3:
            report["root_generation"] = dict(
                radii=list(map(str, radii)), turns=list(map(str, turns)), requested_bits=bits,
                phase_grid=phase_grid, rotation=str(root.rotation), accretivity=str(root.accretivity),
                column=[[str(a), str(b)] for a, b in root.column], column_bits=root.phase_bits,
                phase_bits=phases.bits, phase_internal_bits=phases.working_bits,
            )
        return root, phases

    def fourier_axis(compiler, r, t, direction, bits):
        axis = wrap_axis(compiler, r, t, direction, bits)
        axis_calls.append(dict(true=r, padded=t, direction=direction,
                               requested_bits=bits, generated_layers=len(axis.layers)))
        return axis

    def word(grid, z):
        out = wrap_word(grid, z)
        raw = (z.re, z.im) if hasattr(z, "re") else (Fraction(z), Fraction(0))
        # Keep the actual object alive to preserve identity through lowering.
        origin[id(out)] = (out, raw)
        return out

    def diagnostic_product(self, left, left_bytes, right, right_bytes, output, **kwargs):
        nonlocal machine
        require(machine is None, "unexpected second product")
        machine = self
        def magnitude(tape, count):
            saved = tape.tell()
            data = tape.read(count)
            tape.seek(saved)
            return int.from_bytes(data, "little")
        report.update(left_magnitude=magnitude(left, left_bytes),
                      right_magnitude=magnitude(right, right_bytes),
                      left_negative=kwargs.get("left_negative", False),
                      right_negative=kwargs.get("right_negative", False),
                      sides=self.sides, widths=self.widths, b=self.b, P=self.P, p=self.p,
                      storage_bytes_per_component=self.wide, compiler=axis_calls,
                      pack_spans=[[a, b] for a, b, _ in self.packs],
                      full_one_pass_layer_count=len(self.packs[1][2]))
        require(self.sides == (5, 3) and self.widths == (3, 2) and self.P == 72
                and self.p == 48 and self.b == 8 and self.total == 32, "geometry drift")
        report["direction"] = next(a["direction"] for a in axis_calls if a["true"] == 3)
        return wrap_product(self, left, left_bytes, right, right_bytes, output, **kwargs)

    def crt(tape, total, record, public, route):
        report["before_crt"] = snapshot(tape, machine.wide, total)
        result = wrap_crt(tape, total, record, public, route)
        report["operand_field"] = snapshot(tape, machine.wide, total)
        report["crt_stats"] = dict(machine.tree.stats)
        return result

    def transform(self, tape, opposite=False):
        require(self is machine and not opposite and not finished, "unexpected transform")
        report["opposite"] = opposite
        require(snapshot(tape, self.wide, self.total) == report["operand_field"],
                "post-CRT operand differs at transform entry")
        wrap_transform(self, tape, opposite)
        require(finished and not selected and snapshot(tape, self.wide, self.total) == report["exit"],
                "pack was not completed by return of the forward transform")
        report["transform_returned"] = True
        raise PackComplete

    def diagonal(tape, size, grid, axis_tables):
        nonlocal selected, finished
        active = axis_tables[1]
        if active is machine.entrance[1]:
            require(not selected and not finished, "repeated entry")
            selected = True
            report["entry"] = snapshot(tape, size, machine.total)
            report["h_before_entry"] = machine.h_calls
        if not selected:
            return wrap_diag(tape, size, grid, axis_tables)
        require(axis_tables[0] is machine.idles[0] and all(w == grid.one for w in axis_tables[0]),
                "nonidentity spectator factor")
        values = []
        for pair in active:
            source = origin.get(id(pair))
            rational = (source[1] if source is not None else
                        tuple(Fraction(component, grid.unit) for component in pair))
            # Lowering constructs fixed route signs, ones and halves directly on the grid.
            if source is None:
                require(pair in ((grid.unit, 0), (-grid.unit, 0), (grid.unit // 2, 0)),
                        "unobserved compiler input")
            values.append((tuple(pair), tuple(rational)))
        key = tuple(values)
        if key not in tables:
            tables[key] = len(report["tables"])
            report["tables"].append([dict(word=pair, raw=[str(x) for x in raw])
                                     for pair, raw in key])
        result = wrap_diag(tape, size, grid, axis_tables)
        report["gates"].append(["D", tables[key]])
        report["states"].append(snapshot(tape, size, machine.total))
        if active is machine.exit[1]:
            report["exit"] = report["states"][-1]
            report["h_after_exit"] = machine.h_calls
            finished, selected = True, False
        return result

    def hadamard(tape, widths, positions, public, record, kernel):
        if not selected:
            return wrap_hadamard(tape, widths, positions, public, record, kernel)
        def observed_kernel(src, dst, L, selected_bits, bytes_per_record):
            kernel_calls[(tuple(positions), L, tuple(selected_bits), bytes_per_record)] += 1
            return kernel(src, dst, L, selected_bits, bytes_per_record)
        return wrap_hadamard(tape, widths, positions, public, record, observed_kernel)

    def h(self, tape, start, end, bit):
        result = wrap_h(self, tape, start, end, bit)
        if selected:
            require((start, end) == (1, 2), "unexpected H on another axis")
            report["gates"].append(["H", bit])
            report["states"].append(snapshot(tape, self.wide, self.total))
        return result

    def forbidden(*args, **kwargs):
        raise RuntimeError("unexpected fallback, large native H, or public selector")

    with ExitStack() as hooks:
        for obj, attr, func in (
            (spectral, "polar_root", polar_root),
            (package, "fourier_axis", fourier_axis), (package.Grid, "word", word),
            (product_path.Product, "multiply_tapes", diagnostic_product),
            (product_path, "crt_tape", crt), (multiplier, "tensor_diagonal", diagonal),
            (product_path.Product, "h", h), (product_path, "hadamard_tape", hadamard),
            (product_path, "native_h", forbidden), (multiplier.Product, "transform", transform),
            (multiplier, "stopped_or_retained", forbidden), (multiplier, "multiply_tapes", forbidden),
        ):
            hooks.enter_context(patch.object(obj, attr, func))
        try:
            measure_product.measure(multiplier.Product, 72)
        except PackComplete:
            pass
        require(finished and report.get("transform_returned"), "target pack never returned")

    report["h_kernel_shapes"] = [
        dict(physical_positions=pos, L=L, selected_bits=bits,
             bytes_per_record=size, calls=count)
        for (pos, L, bits, size), count in sorted(kernel_calls.items())
    ]
    # Auxiliary uniquely tagged CRT probe separates physical slots carrying equal values.
    # Exercise only the same CRT, never a second Fourier pack for this observation.
    record = 2 * machine.wide
    tags = BytesIO(b"".join(bytes([i]) * record for i in range(machine.total)))
    wrap_crt(tags, machine.total, record, machine.p, lambda s, d, w: machine.tree.apply(s, d, w))
    routed = tags.getvalue()
    slots = [routed[i:i + record] for i in range(0, len(routed), record)]
    require(all(slot == slot[:1] * record for slot in slots), "torn CRT tag record")
    report["source_slot_at_crt_slot"] = [slot[0] for slot in slots]
    require(sorted(report["source_slot_at_crt_slot"]) == list(range(machine.total)),
            "tagged CRT did not permute the field")
    require(report["operand_field"] == [report["before_crt"][i]
            for i in report["source_slot_at_crt_slot"]], "tagged CRT differs on actual operand")
    return report


def canonical(data):
    return json.dumps(data, sort_keys=True, separators=(",", ":")).encode()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="root of the pinned checkout")
    parser.add_argument("output", type=Path, help="destination .json.gz")
    args = parser.parse_args()
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from provenance import verify
    verify(args.source)
    report = observe(args.source)
    report["payload_sha256"] = hashlib.sha256(canonical(report)).hexdigest()
    with args.output.open("wb") as raw, gzip.GzipFile(filename="", fileobj=raw, mode="wb",
                                                      mtime=0, compresslevel=6) as zipped:
        zipped.write(canonical(report))
    print(f"{PACK}: {len(report['gates'])} gates, {len(report['tables'])} diagonal tables")
    print("payload_sha256", report["payload_sha256"])


if __name__ == "__main__":
    main()
