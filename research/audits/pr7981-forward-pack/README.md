# One completed forward Fourier pack in PR #7981

**Submitted head** `aliromman/integer-mult-bounds@9585c1fa2ed5f2316ac09e04b357cb720badfcb8`  
**Pack ID** `3x5-P72/left-forward/pack-1/physical-axis-1-side-3` (zero-based).

See **RESULTS.md** for the verdict, the minimized witness, the premise map, the control results and the prior-work search. **boundary.tsv** contains all 32 complex entries and exits as pairs of signed integers. Divide each component by $2^{72}$. **pack.json.gz** contains the complete chronological trace, deduplicated coefficient tables, source rationals and CRT tags.

## Run offline

From this directory, run with Python 3.12.

```sh
python -B witness.py pack.json.gz
python -B controls.py pack.json.gz
python -B check.py pack.json.gz
```

The controls command verifies the clean endpoint, doubles the terminal gain on the true rows in a temporary copy, recomputes the affected state and output, refreshes its hash and checks that only the incorrect mean is rejected. The full audit uses exit code 0 for a funded witness, 1 for inconsistent data or an incorrect endpoint, 2 for a certified failure of the spending inequality, and 3 for inconclusive enclosures. The checker needs only the bundled trace, its two Python modules and the standard library.

To retain the altered trace and see both endpoint checks separately, run

```sh
python -B check.py pack.json.gz --endpoint-only
python -B check.py pack.json.gz --negative double-exit-gain.json.gz
python -B check.py double-exit-gain.json.gz --endpoint-only
```

The last command reports successful integrity and gate replay, then returns exit code 1 at physical slot 0.

## Capture again from source

Start from this directory. Check out the submitted head with Git's original LF bytes. Only the package at `research/triangular-native-supplier` is needed.

```sh
git init source
git -C source config core.autocrlf false
git -C source remote add origin https://github.com/aliromman/integer-mult-bounds.git
git -C source sparse-checkout init --cone
git -C source sparse-checkout set research/triangular-native-supplier
git -C source fetch --depth=1 --filter=blob:none origin 9585c1fa2ed5f2316ac09e04b357cb720badfcb8
git -C source checkout --detach FETCH_HEAD
python -I -B capture.py source recaptured.json.gz
python -B check.py recaptured.json.gz
```

The harness checks 70 production blob IDs, including the two-coordinate diagnostic and the Lean sources, and reconstructs the entire package tree ID `2802c261f08c6d42d0f371a1b2adc822e4372d4d`. The production loaders verify their own two nested source archives. For this audit we materialized the pinned package locally. Its outer archive was reconstructed from all 27 checked members and matched the published Git blob `c1a75d8f555f109cb522d3050ce00d75b20efdb8` and SHA-256 `77742c09b088f27e8742c3fad7b7148d3397e4d51b95312106f39a8970165903`.

## What ran

The entry is the existing `measure_product.measure(Product, 72)` with magnitudes 741503 and 466013, marked negative on the left. Its constructor explicitly selects the Fourier geometry. The public small-input multiplier would choose the elementary branch for these operands. The harness wraps the actual routines, restores file cursors after reading, records the constructor and temporary tables, and observes the real stopped Hadamard kernel on the padded field. The first pack runs to supply the entry. Capture ends as soon as the left forward transform returns from its second pack, before truncation to the public grid or the right operand transform. A separate CRT-only tagged permutation distinguishes padded locations containing equal values.

The sorted true extents are `(5,3)` and the padded extents are `(8,4)`. Slot `s+8j` carries spectator coordinate `s` and active coordinate `j`. A digit of ordinal `k` enters the CRT field at `(k mod 5, 2k mod 3)`, since 5 has inverse 2 modulo 3. The initial magnitude has radix-256 digits 127, 80 and 11. Their physical positions are 0, 17 and 10, with values divided by 1024. The pack acts along the ternary coordinate on all eight spectator fibres. There are 15 true addresses and 17 global padding addresses. Even where the spectator is padded, the active coordinate has its own Fourier block. At active coordinate 3 the ideal block copies its input.

The private compiler requests 80 bits. Its three chirp turns are $0,-2/3,-5/3$. The capture records phase grid 6, rotation $-1/12$, a 112-bit root column, 116 phase bits and 296 internal phase bits. It constructs a dyadic circulant root column, synthesizes the padded root, applies its lowered word twice between entrance and exit diagonals, and includes the Gauss-mean compensation on the true exit rows. Each stored complex coefficient uses two exact rationals. Conversion to the 72-bit grid rounds ties to even. Complex data multiplication rounds each component to nearest with ties toward positive infinity. An identity factor preserves its argument exactly. A physical H gate computes unnormalized sums and differences with integer words. Each complex value occupies two signed little-endian fields of 23 bytes each.

## Independent reference

For each spectator fibre with inputs $(z_0,z_1,z_2,z_3)$, the reference returns $(z_0+z_1\omega^k+z_2\omega^{2k})/3$ for $k=0,1,2$, followed by $z_3$, where $\omega=(-1-i\sqrt{3})/2$. This is the negative Fourier character with mean normalization, derived directly from the primitive cube root of unity. `check.py` expands it using rationals and brackets each needed square root between adjacent multiples of $2^{-160}$ using integer square roots. All distance comparisons use the maximum complex modulus across the full padded field. The checker recomputes each transition and coefficient conversion independently of the source program.

The exact gate used to measure each local residual acts with the captured **pre-grid rational coefficient**. This charges factor conversion and data arithmetic to the rounded premise. The synthesis premise compares the unrounded rational gate word with the character formula on this pack's entry. `interval_norm.py` encloses the four-by-four complex matrices on each fibre using dyadic centres and rigorous operator-norm radii. It bounds all contiguous subwords using bounds on every prefix and inverse prefix. This supports the all-inputs Lipschitz statement for this one fixed word. It also encloses the full rational product at the actual pack entry. The compiler's accuracy on other inputs, other packs and other geometries is outside the experiment.

## Credit

**Verifier design credit — Brian Ross.** His Krenn work reconstructs the mathematical object independently and tests semantic corruption after refreshing integrity hashes. His Theta Twist work checks identities with exact coefficient arithmetic. This audit applies those methods to one pack. Ali Romman contributed the completed-pack implementation, integration and Lean arguments. The submitted project's notice credits José Manuel for the rational Fourier compiler. [Krenn](https://github.com/brianross93/Krenn) · [Theta Twist](https://github.com/brianross93/thetatwist) · [Project credit](https://github.com/aliromman/integer-mult-bounds/blob/9585c1fa2ed5f2316ac09e04b357cb720badfcb8/research/triangular-native-supplier/NOTICE.md)
