import time
import converter
import urllib.request

COLORS = {
    "A--": ("gray",   ""),
    "-B-": ("orange", ""),
    "--C": ("cyan",   ""),
    "AB-":("yellow",  ""),
    "A-C":("yellow",  "red"),
    "-BC":("white",  "red"),
    "ABC":("red", ""),
}

def color_map_kegg(map_id, objects, fg, bg):
#	init_kegg()
    # print serv.get_html_of_colored_pathway_by_objects(map_id, objects, fg, bg)
    # print map_id, objects, fg, bg
    # return serv.color_pathway_by_objects(map_id, objects, fg, bg)

    ## first option:
    # retrieve KGML file: http://rest.kegg.jp/get/hsa04666/kgml
    # color it by hand: harder but full freedom!

    ## Alternative:
    # use the existing pathway colorizer:
    # http://www.kegg.jp/kegg-bin/show_pathway?map=hsa00010&multi_query=7167+red%2Cblue%0D00118+pink
    # find the URL of the temp image: <img src="/tmp/...">

    if map_id.startswith("path:"):
        map_id = map_id[5:]

    if len(objects) < 1:
        # get the original image
        return "http://rest.kegg.jp/get/%s/image" % map_id

    data = []
    for i in range(len(objects)):
        o = objects[i]
        ofg = fg[i]
        obg = bg[i]
        data.append("%s+%s,%s" % (o, obg, ofg))

    extra = "%0D".join(data)
    url = "http://www.kegg.jp/kegg-bin/show_pathway?map=%s&multi_query=%s" % (map_id, extra)
    print(url)
    cpage = str( urllib.request.urlopen(url).read() )
    start = cpage.find('<img src="/tmp/')
    start = cpage.find('"', start) + 1
    end = cpage.find('"', start)
    return url,"http://www.kegg.jp" + cpage[start:end]


def prepare_colors(out, d_objects1, d_objects2=set(), d_objects3=set()):
    d_objects = d_objects1.union(d_objects2).union(d_objects3)
    objects = []
    fg = []
    bg = []
    out.write('<table><tr><th>Uniprot</th><th>KEGG</th><th>Symbol</th><th>Value</th><th>MDA D</th><th>MDA sites</th><th>DG75 D</th><th>DG75 sites</th></tr>\n')
    for uid in d_objects:
        kid = converter.handler.to_external(converter.KEGG_IDX, uid)
        if not kid or not kid[0].startswith('hsa:'): continue
        kid = kid[0][4:]
#        b,f = COLORS[d_objects[uid]]
        value = ''
        if uid in d_objects1: value += 'A'
        else: value += '-'
        if uid in d_objects2: value += 'B'
        else: value += '-'
        if uid in d_objects3: value += 'C'
        else: value += '-'
        
        # TODO: show the number of pathways in which the protein is involved (high values means less specificity)
        
        out.write('<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td>' % (uid,kid, converter.handler.to_symbol(uid), value) )
        if uid in node_info:
            for v in node_info[uid]:
                out.write('<td>%s</td>' % v)
        else:
            print('%s is missing' % uid)
            for i in range(4):
                out.write('<td></td>')
        out.write('</tr>\n')
        
        b,f = COLORS[value]
        
        objects.append(kid)
        fg.append(f)
        bg.append(b)
    out.write('</table>\n')
    return objects,fg,bg

# Load pathway names
pathway_names = {}
f = open('SYK_output/mcf7_positif.txt__pathways.tsv')
for line in f:
    data = line.strip('\n').split('\t')
    pathway_names[ data[0] ] = data[4]
f.close()

# Load datasets: list of members for each pathway
filenames = ('SYK_output/mcf7_positif.txt__pathways.tsv','SYK_output/mda231_positif.txt__pathways.tsv','SYK_output/dg75.txt__pathways.tsv',)
datasets = []
for filename in filenames:
    f = open(filename)
    cur_dataset = {}
    for line in f:
        data = line.strip('\n').split('\t')
        pathway = data[0]
        members = data[5].strip()
        if members:
            cur_dataset[pathway] = set(members.split(','))
    f.close()
    datasets.append(cur_dataset)

f = open('SYK_datasets/merged.tsv')
node_info = {}
f.readline()
for line in f:
    data = line.strip('\n').split('\t')
    uid = data[0]
    node_info[uid] = data[2:]
f.close()


# color the selected pathways

# TODO: select the pathways to be colored
out_index = open('kegg_colored/index.html', 'w')
out_index.write("<html>\n<body>\n<ul>\n")
for pathway in pathway_names:
    is_in_all = True
    for d in datasets:
        if pathway not in d:
            is_in_all = False
            break
    if not is_in_all: continue
    
    all_members = [ d[pathway] for d in datasets ]
    out_index.write("<li><a href='%s.html'>%s</a></li>\n" % (pathway,pathway_names[pathway]))
    out = open('kegg_colored/%s.html' % pathway, 'w')
    out.write('<html>\n<head><link href="style.css" rel="stylesheet" type="text/css"></head>\n<body>\n<img src="kegg_colored_legend.png" class="legend" /><br/><img src="%s.png" />\n<pre>' % pathway)
    
    # Replace uids in all_members by gene names
    all_member_names = []
    for mbs in all_members:
        names = [ converter.handler.to_symbol(uid) for uid in mbs ]
        all_member_names.append(names)
    
    out.write('%s\t%s\n' % (pathway, all_member_names))
    objects,fg,bg = prepare_colors(out, *all_members)
    out.write("</pre>\n")
    
    url,img = color_map_kegg(pathway, objects,fg,bg)
    urllib.request.urlretrieve(img, 'kegg_colored/%s.png' % pathway)
    time.sleep(2)

    out.write('<a href="%s">Open in KEGG website</a>\n' % url)
    
    out.write("</ul>\n</body>\n</html>\n")
    out.close()
out_index.write("</body>\n</html>\n")
out_index.close()

