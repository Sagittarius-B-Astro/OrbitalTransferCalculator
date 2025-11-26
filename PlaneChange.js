let pyodide;

async function loadPy() {
    pyodide = await loadPyodide({
        indexURL: "https://cdn.jsdelivr.net/pyodide/v0.26.0/full/"
    });
    await pyodide.loadPackage("numpy");
}

loadPy();

document.getElementById("computeBtn").onclick = async () => {
const result = await pyodide.runPythonAsync(`

    Should contain PlaneChange.py

    `);

document.getElementById("output").textContent = result;
};