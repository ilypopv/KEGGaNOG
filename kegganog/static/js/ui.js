// ==========================================================================
// UI helpers — mode switching, status bar, GitHub stars
// ==========================================================================

function switchMode(mode) {
  currentMode = mode;
  document.getElementById("mode-single").classList.toggle("active", mode === "single");
  document.getElementById("mode-multi").classList.toggle("active",  mode === "multi");
  document.getElementById("left-single").classList.toggle("hidden", mode !== "single");
  document.getElementById("left-multi").classList.toggle("hidden",  mode !== "multi");
  document.getElementById("left-params-single").classList.toggle("hidden", mode !== "single");
  document.getElementById("left-params-multi").classList.toggle("hidden",  mode !== "multi");
  hideResults();
  hideStatus();
}

function showStatus(type, msg) {
  const bar = document.getElementById("status-bar");
  const txt = document.getElementById("status-text");
  const spn = document.getElementById("status-spinner");
  bar.className = `status-bar ${type}`;
  txt.textContent = msg;
  spn.style.display = type === "running" ? "block" : "none";
}

function hideStatus() {
  document.getElementById("status-bar").className = "status-bar";
}

// Live GitHub star count
fetch("https://api.github.com/repos/ilypopv/KEGGaNOG")
  .then(r => r.json())
  .then(d => {
    if (d.stargazers_count !== undefined)
      document.getElementById("gh-stars").textContent = d.stargazers_count;
  })
  .catch(() => {});
