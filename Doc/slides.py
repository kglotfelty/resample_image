
from ciao_contrib.cipt.enhanced_region import *
from chips_contrib.plot_shapes import *


def make_grid( xlen, ylen ):
    
    prop = ChipsRegion()
    prop.depth = 1  # very bottom of plot
    prop.edge.color = "black"
    prop.edge.thickness=1
    prop.fill.style="solid"
    prop.fill.color="white"
    prop.opacity=1
    prop.stem="grid"
    
    for yy in xrange(ylen):
        for xx in xrange(xlen):            
            pixel = rectangle( xx, yy, xx+1,yy+1)
            plot_region(pixel, prop)



def setup_plot( xlen, ylen ):
    clear()
    load_preferences()
    add_window( 240,240, "pixels")
    add_curve([1],[1])
    set_plot_aspect_ratio("1:1")
    set_data_aspect_ratio("1:1")
    delete_curve()
    limits(X_AXIS, -0.5,xlen+0.5)
    limits(Y_AXIS, -0.5,ylen+0.5)
    prop = ChipsPlot()
    prop.bottommargin=0
    prop.topmargin=0
    prop.leftmargin=0
    prop.rightmargin=0
    prop.style="open"
    set_plot(prop)
    hide_axis("ax1")
    hide_axis("ay1")
    
    make_grid(xlen,ylen)




def doit_infile(prop):
    setup_plot( 3,3 )
    pixel = rectangle(1,1,2,2)
    plot_region(pixel, prop)
    print_window("infile.png", "export.clobber=True")


def doit_shifted(prop):
    prop.opacity=1
    setup_plot( 3,3 )
    pixel = rectangle(1.25,1.33,2.25,2.33)
    plot_region(pixel, prop)
    print_window("shifted.png", "export.clobber=True")
    

def doit_area(prop):
    setup_plot(3,3)
    
    p1 = rectangle( 1.25,1.33, 2,2 )
    prop.opacity=p1.area()
    plot_region( p1, prop )

    p2 = rectangle( 1.25, 2, 2, 2.33 )
    prop.opacity=p2.area()
    plot_region( p2, prop )

    p3 = rectangle( 2, 2, 2.25, 2.33 )
    prop.opacity=p3.area()
    plot_region( p3, prop )
    
    p4 = rectangle(2, 1.33, 2.25, 2)
    prop.opacity=p4.area()
    plot_region( p4, prop )
    print_window("area_weight.png", "export.clobber=True")


def doit_reproject(prop):
    setup_plot(3,3)
    
    p1 = rectangle( 1.25,1.33, 2,2 )
    prop.opacity=p1.area()
    p1 = rectangle( 1, 1, 2, 2 )
    plot_region( p1, prop )

    p2 = rectangle( 1.25, 2, 2, 2.33 )
    prop.opacity=p2.area()
    p2 = rectangle( 1, 2, 2, 3 )
    plot_region( p2, prop )

    p3 = rectangle( 2, 2, 2.25, 2.33 )
    prop.opacity=p3.area()
    p3 = rectangle( 2,2, 3,3 )
    plot_region( p3, prop )
    
    p4 = rectangle(2, 1.33, 2.25, 2)
    prop.opacity=p4.area()
    p4 = rectangle(2,1,3,2 )
    plot_region( p4, prop )
    print_window("area_reproject.png", "export.clobber=True")
    

from chips_contrib.lut.hexify import *

def pull_lut( infile ):    
    foo = open(infile,"r").read()
    goo = set(foo.split("\n"))
    moo = [ "".join(map( lambda y: hexify(float(y)), x.split())) for x in goo if len(x) > 0 ]
    return(moo)
    

set_preference( "export.clobber", "True" )
    

def set_prop():
    prop = ChipsRegion()
    prop.depth = 100
    prop.edge.color = "firebrick"
    prop.edge.thickness=1
    prop.fill.style="solid"
    prop.fill.color="firebrick"
    prop.opacity=1
    prop.stem="pixel"
    return prop

prop = set_prop()

lut="/home/kjg/.ds9/LUT/Colo/Colo_Absence.lut"

#doit_infile(prop)
#doit_shifted(prop)
#doit_area(prop)
#doit_reproject(prop)



# ----------------------------------

def doit_alias_in(prop):
    setup_plot(5,5)
    prop.opacity=1

    clrs = pull_lut(lut)
    
    p1 = rectangle(1,1,2,2)
    prop.edge.color = clrs[0]
    prop.fill.color = clrs[0]
    plot_region(p1, prop )

    p2 = rectangle(1,2,2,3)
    prop.edge.color = clrs[1]
    prop.fill.color = clrs[1]
    plot_region(p2, prop )

    p3 = rectangle(2,1,3,2)
    prop.edge.color = clrs[2]
    prop.fill.color = clrs[2]
    plot_region(p3, prop )

    p4 = rectangle(2,2,3,3)
    prop.edge.color = clrs[3]
    prop.fill.color = clrs[3]
    plot_region(p4, prop )

    print_window("map_input.png", "export.clobber=True")
    

def doit_alias_out(prop):
    setup_plot(5,5)
    prop.opacity=1

    clrs = pull_lut(lut)
    
    p1 = rectangle(1,1,2,2)
    prop.edge.color = clrs[0]
    prop.fill.color = clrs[0]
    plot_region(p1, prop )

    p2 = rectangle(1,2,2,3).tweak(dy=1)
    prop.edge.color = clrs[1]
    prop.fill.color = clrs[1]
    plot_region(p2, prop )

    p3 = rectangle(2,1,3,2).tweak(dx=1)
    prop.edge.color = clrs[2]
    prop.fill.color = clrs[2]
    plot_region(p3, prop )

    p4 = rectangle(2,2,3,3).tweak(dx=1,dy=1)
    prop.edge.color = clrs[3]
    prop.fill.color = clrs[3]
    plot_region(p4, prop )

    print_window("map_output.png", "export.clobber=True")
    

def doit_alias_rev(prop):
    setup_plot(5,5)
    prop.opacity=1

    clrs = pull_lut(lut)
    
    prop.edge.color = clrs[0]
    prop.fill.color = clrs[0]
    p1 = rectangle(1,1,2,2)
    plot_region(p1, prop )
    plot_region(p1.tweak(dx=1), prop)
    plot_region(p1.tweak(dy=1), prop)
    plot_region(p1.tweak(dx=1,dy=1), prop)


    prop.edge.color = clrs[1]
    prop.fill.color = clrs[1]
    p2 = rectangle(1,2,2,3).tweak(dy=1)
    plot_region(p2, prop )
    plot_region(p2.tweak(dx=1), prop)
    plot_region(p2.tweak(dy=1), prop)
    plot_region(p2.tweak(dx=1,dy=1), prop)

    prop.edge.color = clrs[2]
    prop.fill.color = clrs[2]
    p3 = rectangle(2,1,3,2).tweak(dx=1)
    plot_region(p3, prop )
    plot_region(p3.tweak(dx=1), prop)
    plot_region(p3.tweak(dy=1), prop)
    plot_region(p3.tweak(dx=1,dy=1), prop)

    prop.edge.color = clrs[3]
    prop.fill.color = clrs[3]
    p4 = rectangle(2,2,3,3).tweak(dx=1,dy=1)
    plot_region(p4, prop )
    plot_region(p4.tweak(dx=1), prop)
    plot_region(p4.tweak(dy=1), prop)
    plot_region(p4.tweak(dx=1,dy=1), prop)

    print_window("map_output_rev.png", "export.clobber=True")
    

def doit_rev_grid(prop, scale=1.0):
    setup_plot(5,5)
    prop.opacity=1
    prop.depth=150
    dd = 0.1
    
    
    c1 = pull_lut("/home/kjg/.ds9/LUT/CB/CB_Paired_12.lut")
    c2 = pull_lut("/home/kjg/.ds9/LUT/CB/CB_Set3_12.lut")
    c1.extend(c2)
    
    ii=0
    for yy in xrange(5):
        yy = (yy/scale)+0.5
        for xx in xrange(5):
            xx=(xx/scale)+0.5
            prop.edge.color = "black"
            prop.fill.color = c1[ii]
            plot_region( rectangle(xx-dd,yy-dd,xx+dd,yy+dd ), prop )
            ii=(ii+1)%len(c1)

    if 2.0 == scale:
        clrs = pull_lut(lut)
        prop.depth=100
        p1 = rectangle(1,1,2,2)
        prop.edge.color = clrs[0]
        prop.fill.color = clrs[0]
        plot_region(p1, prop )

        p2 = rectangle(1,2,2,3)
        prop.edge.color = clrs[1]
        prop.fill.color = clrs[1]
        plot_region(p2, prop )

        p3 = rectangle(2,1,3,2)
        prop.edge.color = clrs[2]
        prop.fill.color = clrs[2]
        plot_region(p3, prop )

        p4 = rectangle(2,2,3,3)
        prop.edge.color = clrs[3]
        prop.fill.color = clrs[3]
        plot_region(p4, prop )

    print_window("map_output_grid_{}.png".format(scale), "export.clobber=True")


#doit_alias_in(prop)
#doit_alias_out(prop)
#doit_alias_rev(prop)
#doit_rev_grid(prop,scale=1.0)
#doit_rev_grid(prop,scale=2.0)




def doit_uniform(prop, dx=0, dy=0, seed=44595):
    setup_plot( 3,3 )
    prop.opacity=1
    pixel = rectangle(1+dx,1+dy,2+dx,2+dy)
    plot_region(pixel, prop)

    prop.depth=150
    dd = 0.07
    npts = 10

    cts = np.zeros(9).reshape([3,3])
    np.random.seed(seed)

    
    prop.edge.color = "black"
    prop.fill.color = "black"
    prop.fill.style = "none"
    for mm in range(npts):
        px = np.random.rand()+1+dx
        py = np.random.rand()+1+dy
        plot_region( rectangle(px-dd,py-dd,px+dd,py+dd ), prop )

        cts[int(py)][int(px)] += 1

    print_window("subpix_{}_{}_{}.png".format(dx,dy,seed), "export.clobber=True")

    prop = set_prop()
    setup_plot(3,3)
    frac = cts/float(np.sum(cts))
    for yy in range( cts.shape[0] ):
        for xx in range(cts.shape[1] ):
            if frac[yy][xx] == 0: continue
            prop.opacity = frac[yy][xx]
            plot_region( rectangle(xx,yy,xx+1,yy+1), prop)

    print_window("brute_{}_{}_{}.png".format(dx,dy,seed), "export.clobber=True")


#prop = set_prop()
#doit_uniform(prop, seed=44595)

#prop = set_prop()
#doit_uniform(prop,0.25,0.33, seed=44595)

    
#prop = set_prop()
#doit_uniform(prop,0.25,0.33, seed=8832)
#prop = set_prop()
#doit_uniform(prop,0.25,0.33, seed=56789)


def doit_prob(prop):
    setup_plot(3,3)
    delete_plot()
    
    p1 = rectangle( 1.25,1.33, 2,2 ).area()

    p2 = rectangle( 1.25, 2, 2, 2.33 ).area()

    p3 = rectangle( 2, 2, 2.25, 2.33 ).area()
    
    p4 = rectangle(2, 1.33, 2.25, 2).area()

    xx = range(4)

    hh = ChipsHistogram()
    hh.fill.color ="firebrick"
    hh.fill.style="solid"

    yy = [p1, 0,0,0]
    hh.fill.opacity=p1
    add_histogram( xx, yy, hh )
    yy = [0, p2, 0,0]
    hh.fill.opacity=p2
    add_histogram( xx, yy, hh )
    yy = [0,0,  p3, 0]
    hh.fill.opacity=p3
    add_histogram( xx, yy, hh )
    yy = [0,0,0, p4]
    hh.fill.opacity=p4
    add_histogram( xx, yy,  hh )



    ax = get_axis("ax1")
    ax.minortick.visible = False    
    set_axis("ax1", ax)
    set_axis("ay1", "offset.perpendicular=30")

    prop = ChipsPlot()
    prop.bottommargin=0.2
    prop.topmargin=0.05
    prop.leftmargin=0.2
    prop.rightmargin=0.05 
    prop.style="open"
    set_plot(prop)
    set_plot_xlabel("Pixel")
    set_plot_ylabel("Probability")

    print_window("pdf.png", "export.clobber=True")

    # ---------------------------------------

    yy = [0, p1, p2, p3, p4]
    print yy
    xx = np.arange(len(yy))-0.5
    xh = xx+1
    cc = np.cumsum(yy)

    add_curve( xx, cc, "symbol.style=none line.thickness=2" )
    delete_histogram("hist1")
    delete_histogram("hist2")
    delete_histogram("hist3")
    delete_histogram("hist4")

    rr = ChipsRegion()
    rr.fill.color = 'firebrick'
    rr.edge.color = 'black'
    rr.edge.style = "none"
    rr.depth=90
    
    cx = np.array([ 0, 1, 1 ])-0.5
    cy = [ cc[0], cc[1], 0 ]
    rr.opacity=p1
    add_region(cx,cy, rr)

    cx = np.array([ 1, 1, 2, 2 ])-0.5
    cy = [ 0, cc[1], cc[2], 0 ]
    rr.opacity=p2
    add_region(cx,cy, rr)
    
    cx = np.array([ 2, 2, 3, 3 ])-0.5
    cy = [ 0, cc[2], cc[3], 0 ]
    rr.opacity=p3
    add_region(cx,cy, rr)

    cx = np.array([ 3,3,4,4 ])-0.5
    cy = [ 0, cc[3], cc[4], 0 ]
    rr.opacity=p4
    add_region(cx,cy, rr)
    

    cx = np.array([ 0, 1, 0 ])-0.5
    cy = [ cc[0], cc[1], cc[1] ]
    rr.opacity=p1
    add_region(cx,cy, rr)

    cx = np.array([ 0, 1, 2, 0 ])-0.5
    cy = [ cc[1], cc[1], cc[2], cc[2] ]
    rr.opacity=p2
    add_region(cx,cy, rr)
    
    cx = np.array([ 0, 2, 3, 0 ])-0.5
    cy = [ cc[2], cc[2], cc[3], cc[3] ]
    rr.opacity=p3
    add_region(cx,cy, rr)

    cx = np.array([ 0, 3, 4, 0 ])-0.5
    cy = [ cc[3], cc[3], cc[4], cc[4] ]
    rr.opacity=p4
    add_region(cx,cy, rr)


    ##for lx,ly in enumerate(cc):
    ##    add_line( -10, ly, lx-0.5, ly, "line.color=gray line.style=shortdash")


    #ccc = [ cc[0], 0, 0, 0]
    #hh.fill.opacity=p1
    #add_histogram(xx,xh,ccc, hh)

    #ccc = [ 0, cc[1], 0, 0]
    #hh.fill.opacity=p2
    #add_histogram(xx,xh,ccc, hh)

    #ccc = [ 0, 0, cc[2], 0]
    #hh.fill.opacity=p3
    #add_histogram(xx,xh,ccc, hh)

    #ccc = [ 0, 0, 0, cc[3]]
    #hh.fill.opacity=p4
    #add_histogram(xx,xh,ccc, hh)


    #set_histogram("dropline=True")
    set_plot_ylabel("Cumulative Probability")

    print_window("cdf.png", "export.clobber=True")

    
    
    

doit_prob(prop)




