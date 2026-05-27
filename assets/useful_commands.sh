# Useful commands
# Convert markdown to HTML with custom font size
pandoc input.md -o output.html -s -H <(echo '<style>html { font-size: 80%; }</style>')