// SPDX-License-Identifier: LGPL-3.0-or-later

const assert = require("node:assert/strict");
const fs = require("node:fs");
const createAUSAXS = require(process.argv[2]);

async function main() {
    const m = await createAUSAXS();
    const allocations = [];
    function alloc(bytes) {
        const pointer = m._malloc(bytes);
        assert.notEqual(pointer, 0);
        allocations.push(pointer);
        return pointer;
    }
    function doubles(values) {
        const pointer = alloc(values.length * 8);
        m.HEAPF64.set(values, pointer / 8);
        return pointer;
    }
    function status() {
        const pointer = alloc(4);
        m.HEAP32[pointer / 4] = -1;
        return pointer;
    }
    function checked(pointer) {
        assert.equal(m.HEAP32[pointer / 4], 0);
    }
    function close(actual, expected) {
        assert.ok(Math.abs(actual - expected) < 1e-3, `${actual} != ${expected}`);
    }

    try {
        const q = doubles([0, 0.1]);
        const x = doubles([0, 10]);
        const y = doubles([0, 0]);
        const z = doubles([0, 0]);
        const w = doubles([1, 1]);
        const intensity = doubles([0, 0]);
        const result = status();
        m._debye_no_ff(q, x, y, z, w, 2, 2, intensity, result);
        checked(result);
        close(m.HEAPF64[intensity / 8], 4);
        close(m.HEAPF64[intensity / 8 + 1], 2 + 2 * Math.sin(1));

        const molecule = m._molecule_from_arrays(x, y, z, w, 2, result);
        checked(result);
        m._molecule_debye_exact_userq(molecule, q, intensity, 2, result);
        checked(result);
        close(m.HEAPF64[intensity / 8], 4);
        close(m.HEAPF64[intensity / 8 + 1], 2 + 2 * Math.sin(1));
        m._molecule_debye_raw_userq(molecule, q, intensity, 2, result);
        checked(result);
        close(m.HEAPF64[intensity / 8], 4);
        close(m.HEAPF64[intensity / 8 + 1], 2 + 2 * Math.sin(1));

        const v = doubles([0, 1]);
        const gx = doubles([0, 0]);
        const gy = doubles([0, 0]);
        const gz = doubles([0, 0]);
        m._molecule_debye_raw_vjp(molecule, q, v, 2, gx, gy, gz, 2, result);
        checked(result);
        close(m.HEAPF64[gx / 8], -m.HEAPF64[gx / 8 + 1]);
        assert.ok(Math.abs(m.HEAPF64[gx / 8]) > 0.01);
        m._deallocate(molecule, result);
        checked(result);

        let output = "";
        const callback = m.addFunction((text, length) => { output += m.UTF8ToString(text, length); }, "vii");
        m._set_output_callback(callback);
        m.ccall("set_setting", null, ["string", "string", "number"], ["threads", "4", result]);
        m._reset_output_callback();
        m.removeFunction(callback);
        checked(result);
        assert.match(output, /only supports one thread/);
        const value = alloc(4);
        const type = alloc(4);
        const setting = m.ccall("get_setting", "number", ["string", "number", "number", "number"], ["threads", value, type, result]);
        checked(result);
        assert.equal(m.UTF8ToString(m.HEAPU32[value / 4]), "1");
        m._deallocate(setting, result);
        checked(result);

        m.FS.writeFile("/input.pdb", fs.readFileSync(process.argv[3]));
        assert.equal(m.ccall("io_is_pdb", "number", ["string", "number"], ["/input.pdb", result]), 1);
        checked(result);
        const pdbMolecule = m.ccall("molecule_from_file", "number", ["string", "number"], ["/input.pdb", result]);
        checked(result);
        m._molecule_debye_raw_userq(pdbMolecule, q, intensity, 2, result);
        checked(result);
        assert.ok(Number.isFinite(m.HEAPF64[intensity / 8]) && m.HEAPF64[intensity / 8] > 0);
        m.ccall("set_setting", null, ["string", "string", "number"], ["exv_model", "average", result]);
        checked(result);
        const modelQ = doubles([0.01, 0.1]);
        m._molecule_debye_userq(pdbMolecule, modelQ, intensity, 2, result);
        checked(result);
        // Native single-threaded build, using a fresh offline residue cache.
        const nativeReference = [6269297.085484935, 3308016.5666994243];
        for (let i = 0; i < nativeReference.length; ++i) {
            const actual = m.HEAPF64[intensity / 8 + i];
            assert.ok(Math.abs(actual - nativeReference[i]) < nativeReference[i] * 1e-3,
                `PDB intensity ${actual} differs from native reference ${nativeReference[i]}`);
        }
        m._deallocate(pdbMolecule, result);
        checked(result);
        m.ccall("get_setting", "number", ["string", "number", "number", "number"], ["not_a_setting", value, type, result]);
        assert.notEqual(m.HEAP32[result / 4], 0);
        const message = alloc(4);
        m._get_last_error_msg(message, result);
        checked(result);
        assert.match(m.UTF8ToString(m.HEAPU32[message / 4]), /Unknown setting/);
        console.log("WebAssembly smoke test passed");
    } finally {
        for (const pointer of allocations) m._free(pointer);
    }
}

main().catch(error => {
    console.error(error);
    process.exitCode = 1;
});
