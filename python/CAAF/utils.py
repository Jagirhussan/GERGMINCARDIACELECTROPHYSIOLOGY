from logging import setLogRecordFactory
import time
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

#Perform causal analysis to get edges
import sklearn
from tigramite import data_processing as pp
from tigramite import plotting as tp
from tigramite.pcmci import PCMCI
from tigramite.independence_tests import ParCorr


def causalAnalysisBasedEdges(timeseries,nodes,g,suffix,alpha=1e-3):
    '''
        Extracts time series samples (every 25 samples over time) at each of the node in nodes.
        A 3x3 stencil is setup around the node and the average time series for this block is used.
        Updates the networkx g with edges that are within the alpha_level of `alpha`
        The tigramite causal graph is saved as `suffix[:-4].png`
    '''
    ts = []
    varnames = []
    n_x = []
    n_y = []
    for i,nd in nodes.items():
        n_x.append(nd[0])
        n_y.append(nd[1])
        blk1 = extract_block(timeseries,nd,(3,3))
        tv = np.zeros(blk1[0,0].shape)
        for i in range(3):
            for j in range(3):
                tv += blk1[i,j].astype(float)
        tv /=9
        ts.append(tv[0:-1:25])
        varnames.append(f"{i}")
    df = np.array(ts).T
    ts = pp.DataFrame(df,var_names=varnames)

    cond_ind_test = ParCorr()
    pcmci = PCMCI(dataframe=ts, cond_ind_test=cond_ind_test)
    results = pcmci.run_pcmci(tau_max=4, pc_alpha=None)
    link_matrix = pcmci.return_significant_links(pq_matrix=results['p_matrix'],
                            val_matrix=results['val_matrix'], alpha_level=alpha)['link_matrix']

   
    tp.plot_graph(
        val_matrix=results['val_matrix'],
        link_matrix=link_matrix,
        var_names=varnames,
        link_colorbar_label='cross-MCI',
        node_colorbar_label='auto-MCI',
        vmin_edges=0.,
        vmax_edges = 0.3,
        edge_ticks=0.05,
        cmap_edges='OrRd',
        vmin_nodes=0,
        vmax_nodes=.5,
        node_ticks=.1,
        cmap_nodes='OrRd',
        node_pos={'x':n_x,'y':n_y}
        ); 
    plt.savefig(f"{suffix[:-4]}.png", format="PNG")

    old_matrix = np.copy(link_matrix)
    slink_matrix = np.zeros(old_matrix.shape, dtype="<U3")
    slink_matrix[:] = ""
    for i, j, tau in zip(*np.where(old_matrix)):
        if tau == 0:
            if old_matrix[j, i, 0] == 0:
                slink_matrix[i, j, 0] = "-->"
                slink_matrix[j, i, 0] = "<--"
            else:
                slink_matrix[i, j, 0] = "o-o"
                slink_matrix[j, i, 0] = "o-o"
        else:
            slink_matrix[i, j, tau] = "-->"

    link_matrix_upper = np.copy(slink_matrix)
    link_matrix_upper[:, :, 0] = np.triu(link_matrix_upper[:, :, 0])

    net = np.any(link_matrix_upper != "", axis=2)

    valm = results["val_matrix"]
    nids = list(nodes.keys())
    
    for i in range(net.shape[0]):
        for j in range(net.shape[1]):
            if i != j:
                if net[i,j]:
                    vm = valm[i,j,:].max()
                    if np.isnan(vm):
                        vm = 1e25 #GERGM does not like nan's, so set it to a large value
                    g.add_edge(nids[i],nids[j],weight=f"{vm}")
            else:
                am = np.argmax(results['val_matrix'][i, j])
                vm = results['val_matrix'][i, j, am]
                if np.isnan(vm):
                    vm = 1e25 #GERGM does not like nan's, so set it to a large value
                pm = results['p_matrix'][i, j, am]
                if np.isnan(pm):
                    pm = 1e25 #GERGM does not like nan's, so set it to a large value
                g.nodes[nids[i]]["vm"] = f"{vm}"
                g.nodes[nids[i]]["pm"] = f"{pm}"

def compute_indices(c, ws, length):
    '''
        Get the block extents
    '''
    # default setting: % operations to accommodate odd/even window sizes
    low, high = c - (ws//2), c + (ws//2) + ws%2 

    # correction for overlap with borders of array
    if low<0:
        low, high = 0, ws
    elif high>length:
        low, high = -ws, None

    return low, high

def extract_block(arr, coords, window_size=(3,3)):
    '''
        Extract a block centered at coords from arr. Block size is window_size
    '''
    # extract array shapes and window sizes into single 
    # variables
    len_r, len_c, _ = arr.shape
    wsr, wsc = window_size

    # extract coords and correct for 0-indexing
    r, c = coords
    r0, c0 = r-1, c-1

    row_low, row_high = compute_indices(r0, wsr, len_r)
    col_low, col_high = compute_indices(c0, wsc, len_c)

    return arr[row_low:row_high, col_low:col_high,:]

def getBlockConnections(connections,index1,blksize=(3,3)):
    '''
        Get the lateral connections at the index, average on a stencil of size blksize
    '''
    blk1 = extract_block(connections,index1,blksize)
    cv = 0  
    for i in range(blksize[0]):
        for j in range(blksize[1]):
            #Count lateral connections
            if blk1[i,j,0]:
                cv +=1
            if blk1[i,j,2]:
                cv +=1
    return cv

def generateCausalNetwork(tissue,filename,rc,cc,alpha=1e-3):
    '''
    Generate a sub network that uniformly samples the actual tissue
    The cells along a row are connected, cells along the first column are connected
    First column is the first column of the tissue
    '''
    selRows = list(range(0,tissue.row_size,tissue.row_size//(rc+2)))[1:rc+1]
    selCols = range(0,tissue.col_size,tissue.col_size//cc)
    row_size = len(selRows)
    col_size = len(selCols)
    connections = tissue.connections
    defects = tissue.defects
    state_timeseries = tissue.state_timeseries
    cell_refractory_period = tissue.cell_refractory_period
    conduction_velocity = tissue.conduction_velocity    
    #Create the graph
    g = nx.DiGraph()    
    #Add nodes
    pos =dict()
    for i,r in enumerate(selRows):
        for j,c in enumerate(selCols):
            nid = i*col_size + j + 1
            cv = 0.0
            if conduction_velocity[r,c,1]!=0:
                cv = conduction_velocity[r,c,0]/conduction_velocity[r,c,1]/c
            pos[nid] = (r,c)
            cor = getBlockConnections(connections,(r,c))
            g.add_node(nid,x=f"{c}",y=f"{r}",rft=f"{cell_refractory_period[r,c]}",defect=f"{defects[r,c]}",cv=f"{cv}",con=f"{cor}")


    causalAnalysisBasedEdges(state_timeseries,pos,g,filename+"_trig",alpha)
    nx.write_pajek(g,filename)
