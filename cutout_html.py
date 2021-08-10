import abell_cluster_module as ab
import webbrowser

fn = ab.work_dir+'pics/A754_cutout/A754.html'

f = open(fn,'w')

message = """<html>
<head>
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>
* {
  box-sizing: border-box;
}

.row {
  display: flex;
}

/* Create three equal columns that sits next to each other */
.column {
  flex: 25%;
  padding: 5px;
}

</style>
</head>
<body>

<h2>Abell 754 member (WINGS-SPE) galaxies</h2>
<p>Each row is for a galaxy in u, g, r, and color.</p>
<p>Color image is a composite of u, g, and r in blue, green, and red, respectively.</p>
"""

for i in range(1, 108):
    div = f"""
    <div class="row">
        <div class="column">
        <img src="pics/{i}_u.png" alt="u" style="width:100%">
      </div>
      <div class="column">
        <img src="pics/{i}_g.png" alt="g" style="width:100%">
      </div>
      <div class="column">
        <img src="pics/{i}_r.png" alt="r" style="width:100%">
      </div>
      <div class="column">
        <img src="pics/{i}_rgb.png" alt="rgb" style="width:100%">
      </div>
    </div>
    """
    message+=div

message+="""
</body>
</html>"""

f.write(message)
f.close()

filename = 'file://'+fn
webbrowser.open_new_tab(filename)