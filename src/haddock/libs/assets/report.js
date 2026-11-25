// Script of the analysis report.html, inlined by haddock.libs.libplots.
//
// Sorts the tables on a click on a header, again to reverse; for tables
// with one column per record (row headers), the columns are reordered.
// Opens the models in the Mol* viewer, loaded on first use from the URL
// in the data-molstar attribute of the viewer dialog; models may be gzipped.

function sortValue(cell) {
  const value = cell.dataset.value;
  const x = parseFloat(value);
  return isNaN(x) ? value : x;
}

document.querySelectorAll('th.sortable').forEach(th => {
  th.addEventListener('click', () => {
    const table = th.closest('table');
    const asc = th.dataset.dir !== 'asc';
    table.querySelectorAll('th.sortable').forEach(h => delete h.dataset.dir);
    th.dataset.dir = asc ? 'asc' : 'desc';
    // Cells without value (NA) go last.
    const compare = (a, b) => {
      const x = sortValue(a), y = sortValue(b);
      if (x === undefined || y === undefined) {
        return (x === undefined) - (y === undefined);
      }
      return (asc ? 1 : -1) * (x < y ? -1 : x > y ? 1 : 0);
    };
    if (th.scope === 'row') {
      const cells = [...th.parentElement.cells].slice(1);
      const order = cells.map((_, i) => i)
        .sort((i, j) => compare(cells[i], cells[j]));
      for (const row of table.rows) {
        const rowCells = [...row.cells].slice(1);
        order.forEach(i => row.appendChild(rowCells[i]));
      }
    } else {
      const body = table.tBodies[0];
      [...body.rows]
        .sort((a, b) => compare(a.cells[th.cellIndex], b.cells[th.cellIndex]))
        .forEach(row => body.appendChild(row));
    }
  });
});

async function createViewer(base) {
  const css = document.createElement('link');
  css.rel = 'stylesheet';
  css.href = `${base}/molstar.css`;
  document.head.appendChild(css);
  await new Promise((resolve, reject) => {
    const script = document.createElement('script');
    script.src = `${base}/molstar.js`;
    script.onload = resolve;
    script.onerror = () => reject(new Error(`cannot load ${script.src}`));
    document.head.appendChild(script);
  });
  return molstar.Viewer.create('viewer', {
    layoutIsExpanded: false,
    layoutShowControls: true,
    layoutShowLog: false,
    collapseLeftPanel: true,
    // no list of states from the Mol* server
    layoutShowRemoteState: false,
    viewportShowExpand: true,
  });
}

let viewer;
async function view(url) {
  const dialog = document.getElementById('viewer-dialog');
  const link = document.getElementById('viewer-file');
  const error = document.getElementById('viewer-error');
  link.href = url;
  link.download = link.textContent = url.split('/').pop();
  error.textContent = '';
  dialog.showModal();
  try {
    viewer ??= createViewer(dialog.dataset.molstar);
    const molViewer = await viewer;
    await molViewer.plugin.clear();
    const response = await fetch(url);
    if (!response.ok) {
      throw new Error(`${response.status} ${response.statusText}`);
    }
    let data = new Uint8Array(await response.arrayBuffer());
    // A server may already have decompressed a .gz file, so look at
    // the gzip magic number rather than at the file name.
    if (data[0] === 0x1f && data[1] === 0x8b) {
      const stream = new Blob([data]).stream()
        .pipeThrough(new DecompressionStream('gzip'));
      data = new Uint8Array(await new Response(stream).arrayBuffer());
    }
    await molViewer.loadStructureFromData(
      new TextDecoder().decode(data), 'pdb');
  } catch (e) {
    error.textContent = ` Could not load the structure: ${e.message}`;
  }
}

document.querySelectorAll('button.view').forEach(button => {
  button.addEventListener('click', () => view(button.dataset.url));
});
