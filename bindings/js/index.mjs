// @astrologer/ephemera -- entry point.  Re-exports the wrapper and adds
// createEphemera(), which finds the packaged .wasm modules and data files
// next to this file (Node: read from disk; browsers and bundlers: fetched).
//
//   import { createEphemera, Body } from '@astrologer/ephemera';
//   const eph = await createEphemera();                      // 1" data built in
//   const eph = await createEphemera({ data: ['stars'] });   // + fixed stars
//   const eph = await createEphemera({ module: 'core', data: ['01as', 'chiron'] });

import { fromWasm } from './ephemeris.mjs';

export * from './ephemeris.mjs';

const MODULES = {
  '1as': 'ephemeris-1as.wasm',               // 1" series built in
  '1as-chiron': 'ephemeris-1as-chiron.wasm', // ... with Chiron
  core: 'ephemeris.wasm',                    // no data: load a file
  lite: 'ephemeris-lite.wasm',               // core without eclipses and stars
};

async function readAsset(path) {
  const url = new URL(path, import.meta.url);
  if (url.protocol === 'file:') {
    const { readFile } = await import('node:fs/promises');
    return new Uint8Array(await readFile(url));
  }
  const response = await fetch(url);
  if (!response.ok) throw new Error(`cannot fetch ${url}: ${response.status}`);
  return new Uint8Array(await response.arrayBuffer());
}

export async function createEphemera({ module = '1as', data = [] } = {}) {
  const file = MODULES[module];
  if (!file) throw new Error(`unknown module "${module}"`);
  const eph = await fromWasm(await readAsset(`./wasm/${file}`));
  for (const d of data) {
    const bytes = typeof d === 'string' ? await readAsset(`./data/eph-${d}.bin`) : d;
    if (!eph.loadData(bytes)) throw new Error(`cannot load data ${typeof d === 'string' ? d : ''}`);
  }
  return eph;
}

/** @deprecated Renamed createEphemera (0.4.0). */
export const load = createEphemera;
