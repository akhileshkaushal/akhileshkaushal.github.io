document.addEventListener('DOMContentLoaded', () => {
  const params = new URLSearchParams(window.location.search);
  const id = params.get('id');

  if (!id) {
    document.body.innerHTML = "<p class='text-danger'>Error: No pipeline ID provided.</p>";
    return;
  }

  fetch('pipelines.json')
    .then(response => response.json())
    .then(data => {
      const pipeline = data.find(p => p.id === id);
      if (!pipeline) {
        document.body.innerHTML = `<p class='text-danger'>Error: No pipeline found for ID '${id}'</p>`;
        return;
      }

      document.title = `${pipeline.name} | Pipeline Details`;
      document.getElementById('pipeline-name').textContent = pipeline.name;

      document.getElementById('pipeline-description').innerHTML = `<p>${pipeline.description}</p>`;
      document.getElementById('pipeline-tags').innerHTML = (pipeline.tags || []).map(tag =>
        `<span class="tag">${tag}</span>`).join('');

      document.getElementById('pipeline-steps').innerHTML = (pipeline.details.steps || []).map(
        step => `<li>${step}</li>`
      ).join('');

      document.getElementById('pipeline-repo').href = pipeline.details.repo;
      document.getElementById('pipeline-version').textContent = `Version: ${pipeline.version}`;
      document.getElementById('pipeline-usage').textContent = pipeline.details.usage;
    })
    .catch(err => {
      console.error(err);
      document.body.innerHTML = "<p class='text-danger'>Failed to load pipeline data.</p>";
    });
});
