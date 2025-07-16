document.addEventListener("DOMContentLoaded", function () {
  fetch("partials/header.html")
    .then(response => response.text())
    .then(data => {
      const header = document.createElement("div");
      header.innerHTML = data;
      document.body.insertBefore(header, document.body.firstChild);
    })
    .catch(error => {
      console.error("Failed to load header:", error);
    });
});
