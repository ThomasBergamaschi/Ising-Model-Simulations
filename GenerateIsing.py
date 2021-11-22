import numpy as np
from scipy import stats


def lattice(d,N):
    '''
    Generates a hypercubic lattice of dimension d with side length N
    '''    
    return 2*np.random.randint(2,size=np.int64((N*np.ones(d))))-1

def energy_change(grid,K,N,d,selected_spin):
    '''
    Computes the change in energy from flipping the selected_spin
    in the lattice grid
    '''
    energy_diff = 0
    for i in range(d):
        neighbor_vec = np.zeros(d)
        neighbor_vec[i] = 1
        energy_diff += 2*K*(grid[tuple(selected_spin)])*grid[tuple(np.int64((selected_spin+neighbor_vec)%N))]
        energy_diff += 2*K*(grid[tuple(selected_spin)])*grid[tuple(np.int64((selected_spin-neighbor_vec)%N))]
        
    return energy_diff


def metropolis(grid,N,d,M,K,mag=False):
    '''
    Applies the metropolis algorithm to the d-dimensional lattice grid of 
    size N  with M repetitions and hamiltonian parameter K.
    '''
    mags = []
    
    for i in range(M):
        selected_spin = np.random.randint(0,N,size=d)
        
        if i>M/2:
            mags.append(np.abs(np.mean(grid)))

        energy_diff = energy_change(grid,K,N,d,selected_spin)

        if energy_diff<0:
            grid[tuple(selected_spin)] = -grid[tuple(selected_spin)]
        else:
            if np.random.random()<np.exp(-energy_diff):
                grid[tuple(selected_spin)] = -grid[tuple(selected_spin)]


        
    if mag:
        return grid,mags
    return grid




def sweep(N,d,M,K_low,dK,steps,iters):
    '''
    Obtains the expectation value of the magnetization
    as a function of K from K_low to K_low + dK*steps
    by simulating iters time each value of K to get
    statistics
    '''
    K = K_low
    K_values = []
    magnetizations = []
    stds = []
    for i in range(steps):
        mags = []
        K_values.append(K)
        for j in range(iters):
            
            grid = lattice(d,N)

            grid,mag = metropolis(grid,N,d,M,K,True)

            mags.append(np.mean(mag))#average over time
            
        K = K + dK
        magnetizations.append(np.mean(mags))#average over trials
        stds.append(np.std(mags))

    return K_values,np.abs(magnetizations),stds

def N_behaviour(N_list,d,K_low,dK,steps,iters):
    """
    Obtains the behaviour of the magnetization as a function
    of the number of spins N in N_list
    """
    magnetizations = []
    magnetizations_std = []
    K_values = []
    
    
    
    for N in N_list:
        M = 2*(N**d)*1000
        k_vals,mags,stds = sweep(N,d,M,K_low,dK,steps,iters)
        
        K_values.append(k_vals)
        magnetizations.append(mags)
        magnetizations_std.append(stds)
        
    return N_list, K_values, magnetizations,magnetizations_std
    
def wolff(grid,N,d,M,K,mag=False):
    '''
    Applies the wolff algorithm to the d-dimensional lattice grid of 
    size N  with M repetitions and hamiltonian parameter K.
    '''
    corrs = []

    middle_coord = np.int64(int(N/2)*np.ones(d))
    middle = tuple(np.int64((N/2)*np.ones(d)))
    
    for i in range(M):
        spin_stack = []

        selected_bonds = []

        cluster = []

        selected_spin = np.random.randint(0,N,size=d)

        cluster.append(tuple(selected_spin))
        
        j = 0

        while len(spin_stack)!=0 or j==0:

            update_cluster(N,d,K,grid,selected_spin,cluster,selected_bonds,spin_stack)

            if len(spin_stack)==0:
                break

            selected_spin = spin_stack.pop()
            j = j+1
        
        if i>3*M/4:
            
            row_corr = []
            nbr_vec = np.ones(d)

            for i in range(int(N/2-1)):
                row_corr.append(grid[middle]*grid[tuple(np.int64(middle_coord+i*nbr_vec))])

            corrs.append(row_corr)


        flip_cluster(grid,cluster)

    if mag:
        return grid,np.mean(np.array(corrs),axis=0)
    return grid
    
def update_cluster(N,d,K,grid,selected_spin,cluster,selected_bonds,spin_stack):
    
    
    for i in range(d):
        neighbor_vec = np.zeros(d)
        neighbor_vec[i] = 1

        c1 = ((grid[tuple(np.int64((selected_spin+neighbor_vec)%N))]*grid[tuple(selected_spin)] == 1))

        c2 = (((tuple(np.int64((selected_spin+neighbor_vec)%N)),tuple(selected_spin)) not in selected_bonds))

        c3 = (((tuple(selected_spin),tuple(np.int64((selected_spin+neighbor_vec)%N))) not in selected_bonds))


        if c1 and c2 and c3:
            

            c4 = (tuple(np.int64((selected_spin+neighbor_vec)%N))) not in cluster

            if (np.random.rand()<1-np.exp(-2*K)) and c4:
                cluster.append(tuple(np.int64((selected_spin+neighbor_vec)%N)))
                spin_stack.append(np.int64((selected_spin+neighbor_vec)%N))


        c1 = ((grid[tuple(np.int64((selected_spin-neighbor_vec)%N))]*grid[tuple(selected_spin)] == 1))

        c2 = (((tuple(np.int64((selected_spin-neighbor_vec)%N)),tuple(selected_spin)) not in selected_bonds))

        c3 = (((tuple(selected_spin),tuple(np.int64((selected_spin-neighbor_vec)%N))) not in selected_bonds))

        if c1 and c2 and c3:

            c4 = (tuple(np.int64((selected_spin-+neighbor_vec)%N))) not in cluster
            
            if np.random.rand()<1-np.exp(-2*K) and c4:
                cluster.append(tuple(np.int64((selected_spin-neighbor_vec)%N)))
                spin_stack.append(tuple(np.int64((selected_spin-neighbor_vec)%N)))

        selected_bonds.append((tuple(np.int64((selected_spin+neighbor_vec)%N)),tuple(selected_spin)))
        selected_bonds.append((tuple(np.int64((selected_spin-neighbor_vec)%N)),tuple(selected_spin)))





def flip_cluster(grid,cluster):
    #print (len(cluster))
    flipped = []
    for spin in cluster:        
        if not(tuple(spin) in flipped):
            grid[spin] = -grid[spin]
            flipped.append(spin)





def correlations(N,d,M,K,iters):
    '''
    Obtains the expectation value of the magnetization
    as a function of K from K_low to K_low + dK*steps
    by simulating iters time each value of K to get
    statistics
    '''

    stds = []
    corrs = []
    correlation_list = []
    for j in range(iters):
        
        grid = lattice(d,N)

        grid,corr = wolff(grid,N,d,M,K,True)#get correlation

        #print (corr)

        corrs.append(corr)#average over time
        
    correlation_list.append(np.mean(np.array(corrs),axis=0))#average over trials

    stds.append(np.std(np.array(corrs),axis=0))

    return correlation_list,stds




#Potts Model

def potts_grid(q,N):
    q_list = list(range(1,q+1))
    return np.random.choice(q_list,size=(N,N))

def potts(grid,K,q,N,M,mags = False):
    '''
    Applies the Wolff algorithm to the potts Model
    with q-valued spins
    '''
    magnetizations = []
    cond = True
    index = 0

    count = 0
    for i in range(M):
        selected_spin = (np.random.randint(0,N),np.random.randint(0,N))
        bonds = {}
        cluster = set()
        stack = set()
        j = 0

        q_list = list(range(1,q+1))

        while len(stack)>0 or j==0:
            cluster.add(selected_spin)
            j = j+1
            nbr1 = (selected_spin[0],(selected_spin[1]+1)%N)
            nbr2 = (selected_spin[0],(selected_spin[1]-1)%N)

            if grid[selected_spin]==grid[nbr1] and np.random.rand()<1-np.exp(-4*K):
                if selected_spin in bonds:
                    if nbr1 not in bonds[selected_spin]:
                        cluster.add(nbr1)
                        stack.add(nbr1)
                else:
                    cluster.add(nbr1)
                    stack.add(nbr1)

            if grid[selected_spin]==grid[nbr2] and np.random.rand()<1-np.exp(-4*K):
                if selected_spin in bonds:
                    if nbr2 not in bonds[selected_spin]:
                        cluster.add(nbr2)
                        stack.add(nbr2)
                else:
                    cluster.add(nbr2)
                    stack.add(nbr2)

            update_bonds(bonds,selected_spin,nbr1,nbr2)
            
            nbr1 = ((selected_spin[0]+1)%N,selected_spin[1])
            nbr2 = ((selected_spin[0]-1)%N,selected_spin[1])

            if grid[selected_spin]==grid[nbr1] and np.random.rand()<1-np.exp(-4*K):
                if selected_spin in bonds:
                    if nbr1 not in bonds[selected_spin]:
                        cluster.add(nbr1)
                        stack.add(nbr1)
                else:
                    cluster.add(nbr1)
                    stack.add(nbr1)

            if grid[selected_spin]==grid[nbr2] and np.random.rand()<1-np.exp(-4*K):
                if selected_spin in bonds:
                    if nbr2 not in bonds[selected_spin]:
                        cluster.add(nbr2)
                        stack.add(nbr2)
                else:
                    cluster.add(nbr2)
                    stack.add(nbr2)

            update_bonds(bonds,selected_spin,nbr1,nbr2)
            if len(stack)==0:
                break

            selected_spin = stack.pop()

        random_cluster(q_list,grid,cluster)
        if len(cluster)>0.99*N*N:
            count = count + 1
        if count>100:
            break
        if count ==1:
            index = i
        
        magnetizations.append(((stats.mode(grid,axis=None)).count)[0])
    
    if mags:
        return grid,np.mean(magnetizations[max(8*len(magnetizations)//9,index):])
    return grid
            
def update_bonds(bonds,selected_spin,nbr1,nbr2):
    if selected_spin in bonds:
        bonds[selected_spin].add(nbr1)
        bonds[selected_spin].add(nbr2)
    else:
        bonds[selected_spin] = set()
        bonds[selected_spin].add(nbr1)
        bonds[selected_spin].add(nbr2)

    if nbr1 in bonds:
        bonds[nbr1].add(selected_spin)
    else:
        bonds[nbr1] = set()
        bonds[nbr1].add(selected_spin)

    if nbr2 in bonds:
        bonds[nbr2].add(selected_spin)
    else:
        bonds[nbr2] = set()
        bonds[nbr2].add(selected_spin)

def random_cluster(q_list,grid,cluster):
    val = np.random.choice(q_list)
    for spin in cluster:
        grid[spin] = val

def potts_magnetization(grid,Klow,dK,steps,q,N,M,iters):
    K = Klow
    magnetizations = []
    stds = []
    K_list = []

    for i in range(steps):
        K_list.append(K)
        mags = []


        for _ in range(iters):

            
            grid = potts_grid(q,N)
            grid, mag = potts(grid,K,q,N,M,True)
            mags.append(mag)


        K = K + dK

        magnetizations.append(np.mean(mags)/(N*N))
        stds.append((np.std(mags))/(N*N))

    return K_list,magnetizations,stds


def ising_3d(grid,K,N,M,mags = False):
    '''
    Applies the Wolff algorithm to a 3d ising model
    '''
    magnetizations = []
    cond = True
    index = 0

    count = 0
    for i in range(M):
        selected_spin = (np.random.randint(0,N),np.random.randint(0,N),np.random.randint(0,N))
        bonds = {}
        cluster = set()
        stack = set()
        j = 0

        while len(stack)>0 or j==0:
            cluster.add(selected_spin)
            j = j+1
            nbr1 = (selected_spin[0],(selected_spin[1]+1)%N,selected_spin[2])
            nbr2 = (selected_spin[0],(selected_spin[1]-1)%N,selected_spin[2])

            if grid[selected_spin]==grid[nbr1] and np.random.rand()<1-np.exp(-2*K):
                if selected_spin in bonds:
                    if nbr1 not in bonds[selected_spin]:
                        cluster.add(nbr1)
                        stack.add(nbr1)
                else:
                    cluster.add(nbr1)
                    stack.add(nbr1)

            if grid[selected_spin]==grid[nbr2] and np.random.rand()<1-np.exp(-2*K):
                if selected_spin in bonds:
                    if nbr2 not in bonds[selected_spin]:
                        cluster.add(nbr2)
                        stack.add(nbr2)
                else:
                    cluster.add(nbr2)
                    stack.add(nbr2)

            update_bonds(bonds,selected_spin,nbr1,nbr2)
            
            nbr1 = ((selected_spin[0]+1)%N,selected_spin[1],selected_spin[2])
            nbr2 = ((selected_spin[0]-1)%N,selected_spin[1],selected_spin[2])

            if grid[selected_spin]==grid[nbr1] and np.random.rand()<1-np.exp(-2*K):
                if selected_spin in bonds:
                    if nbr1 not in bonds[selected_spin]:
                        cluster.add(nbr1)
                        stack.add(nbr1)
                else:
                    cluster.add(nbr1)
                    stack.add(nbr1)

            if grid[selected_spin]==grid[nbr2] and np.random.rand()<1-np.exp(-2*K):
                if selected_spin in bonds:
                    if nbr2 not in bonds[selected_spin]:
                        cluster.add(nbr2)
                        stack.add(nbr2)
                else:
                    cluster.add(nbr2)
                    stack.add(nbr2)
                    
            update_bonds(bonds,selected_spin,nbr1,nbr2)
                    
            #3d Change
                    
            nbr1 = ((selected_spin[0])%N,selected_spin[1],(selected_spin[2]+1)%N)
            nbr2 = ((selected_spin[0])%N,selected_spin[1],(selected_spin[2]-1)%N)

            if grid[selected_spin]==grid[nbr1] and np.random.rand()<1-np.exp(-2*K):
                if selected_spin in bonds:
                    if nbr1 not in bonds[selected_spin]:
                        cluster.add(nbr1)
                        stack.add(nbr1)
                else:
                    cluster.add(nbr1)
                    stack.add(nbr1)

            if grid[selected_spin]==grid[nbr2] and np.random.rand()<1-np.exp(-2*K):
                if selected_spin in bonds:
                    if nbr2 not in bonds[selected_spin]:
                        cluster.add(nbr2)
                        stack.add(nbr2)
                else:
                    cluster.add(nbr2)
                    stack.add(nbr2)

            update_bonds(bonds,selected_spin,nbr1,nbr2)
            if len(stack)==0:
                break

            selected_spin = stack.pop()

        random_cluster_3d(grid,cluster)
        if len(cluster)>0.99*N*N*N:
            count = count + 1
        if count>100:
            break
        if count ==1:
            index = i
        
        magnetizations.append(np.mean(grid))
    
    if mags:
        magnetizations = np.abs(magnetizations)
        return grid,np.mean(magnetizations[max(8*len(magnetizations)//9,index):])
    return grid
    
def random_cluster_3d(grid,cluster):
    for spin in cluster:
        grid[spin] = -grid[spin]
                              
def ising3_magnetization(grid,Klow,dK,steps,N,M,iters):
    K = Klow
    magnetizations = []
    stds = []
    K_list = []

    for i in range(steps):
        K_list.append(K)
        mags = []


        for _ in range(iters):

            
            grid = lattice(3,N)
            grid, mag = ising_3d(grid,K,N,M,True)
            mags.append(mag)


        K = K + dK

        magnetizations.append(np.mean(mags))
        stds.append((np.std(mags)))

    return grid,K_list,magnetizations,stds

def ising2_correlations(grid,K,N,M,iters):
    correlations = []
    stds = []

    for _ in range(iters):

        grid = lattice(2,N)
        grid, corr = ising_2d(grid,K,N,M,True)
        correlations.append(corr)

        stds.append(corr)

    return np.mean(correlations,axis=0),np.std(stds,axis=0)

def ising_2d(grid,K,N,M,corrs = False):
    '''
    Applies the Wolff algorithm to a 2d ising model
    '''
    correlations = []
    cond = True
    index = 0

    count = 0
    for i in range(M):
        selected_spin = (np.random.randint(0,N),np.random.randint(0,N))
        bonds = {}
        cluster = set()
        stack = set()
        j = 0

        while len(stack)>0 or j==0:
            cluster.add(selected_spin)
            j = j+1
            nbr1 = (selected_spin[0],(selected_spin[1]+1)%N)
            nbr2 = (selected_spin[0],(selected_spin[1]-1)%N)

            if grid[selected_spin]==grid[nbr1] and np.random.rand()<1-np.exp(-2*K):
                if selected_spin in bonds:
                    if nbr1 not in bonds[selected_spin]:
                        cluster.add(nbr1)
                        stack.add(nbr1)
                else:
                    cluster.add(nbr1)
                    stack.add(nbr1)

            if grid[selected_spin]==grid[nbr2] and np.random.rand()<1-np.exp(-2*K):
                if selected_spin in bonds:
                    if nbr2 not in bonds[selected_spin]:
                        cluster.add(nbr2)
                        stack.add(nbr2)
                else:
                    cluster.add(nbr2)
                    stack.add(nbr2)

            update_bonds(bonds,selected_spin,nbr1,nbr2)
            
            nbr1 = ((selected_spin[0]+1)%N,selected_spin[1])
            nbr2 = ((selected_spin[0]-1)%N,selected_spin[1])

            if grid[selected_spin]==grid[nbr1] and np.random.rand()<1-np.exp(-2*K):
                if selected_spin in bonds:
                    if nbr1 not in bonds[selected_spin]:
                        cluster.add(nbr1)
                        stack.add(nbr1)
                else:
                    cluster.add(nbr1)
                    stack.add(nbr1)

            if grid[selected_spin]==grid[nbr2] and np.random.rand()<1-np.exp(-2*K):
                if selected_spin in bonds:
                    if nbr2 not in bonds[selected_spin]:
                        cluster.add(nbr2)
                        stack.add(nbr2)
                else:
                    cluster.add(nbr2)
                    stack.add(nbr2)
                    
            update_bonds(bonds,selected_spin,nbr1,nbr2)
                    
            if len(stack)==0:
                break

            selected_spin = stack.pop()

        random_cluster_3d(grid,cluster)
        if len(cluster)>0.99*N*N:
            count = count + 1
        if count>100:
            break
        if count ==1:
            index = i
        
        row_corr = []
        nbr_vec = np.array([1,0])
        middle_coord = np.ones(2)*int(N/2)
        middle = (N//2,N//2)

        for i in range(int(N/2-1)):
            row_corr.append(grid[middle]*grid[tuple(np.int64(middle_coord+i*nbr_vec))])
           
        correlations.append(row_corr)
    
    if corrs:
        return grid,np.mean(correlations[max(5*len(correlations)//9,index):],axis=0)
    return grid


def ising2_magnetization(grid,Klow,dK,steps,N,M,iters):
    K = Klow
    magnetizations = []
    stds = []
    K_list = []

    for i in range(steps):
        K_list.append(K)
        mags = []


        for _ in range(iters):

            
            grid = lattice(2,N)
            grid, mag = ising_2d_magnetizations(grid,K,N,M,True)
            mags.append(mag)


        K = K + dK

        magnetizations.append(np.mean(mags))
        stds.append((np.std(mags)))

    return K_list,magnetizations,stds

def ising_2d_magnetizations(grid,K,N,M,mags = False):
    '''
    Applies the Wolff algorithm to a 2d ising model
    '''
    magnetizations = []
    cond = True
    index = 0

    count = 0
    for i in range(M):
        selected_spin = (np.random.randint(0,N),np.random.randint(0,N))
        bonds = {}
        cluster = set()
        stack = set()
        j = 0

        while len(stack)>0 or j==0:
            cluster.add(selected_spin)
            j = j+1
            nbr1 = (selected_spin[0],(selected_spin[1]+1)%N)
            nbr2 = (selected_spin[0],(selected_spin[1]-1)%N)

            if grid[selected_spin]==grid[nbr1] and np.random.rand()<1-np.exp(-2*K):
                if selected_spin in bonds:
                    if nbr1 not in bonds[selected_spin]:
                        cluster.add(nbr1)
                        stack.add(nbr1)
                else:
                    cluster.add(nbr1)
                    stack.add(nbr1)

            if grid[selected_spin]==grid[nbr2] and np.random.rand()<1-np.exp(-2*K):
                if selected_spin in bonds:
                    if nbr2 not in bonds[selected_spin]:
                        cluster.add(nbr2)
                        stack.add(nbr2)
                else:
                    cluster.add(nbr2)
                    stack.add(nbr2)

            update_bonds(bonds,selected_spin,nbr1,nbr2)
            
            nbr1 = ((selected_spin[0]+1)%N,selected_spin[1])
            nbr2 = ((selected_spin[0]-1)%N,selected_spin[1])

            if grid[selected_spin]==grid[nbr1] and np.random.rand()<1-np.exp(-2*K):
                if selected_spin in bonds:
                    if nbr1 not in bonds[selected_spin]:
                        cluster.add(nbr1)
                        stack.add(nbr1)
                else:
                    cluster.add(nbr1)
                    stack.add(nbr1)

            if grid[selected_spin]==grid[nbr2] and np.random.rand()<1-np.exp(-2*K):
                if selected_spin in bonds:
                    if nbr2 not in bonds[selected_spin]:
                        cluster.add(nbr2)
                        stack.add(nbr2)
                else:
                    cluster.add(nbr2)
                    stack.add(nbr2)
                    
            update_bonds(bonds,selected_spin,nbr1,nbr2)
                    
            if len(stack)==0:
                break

            selected_spin = stack.pop()

        random_cluster_3d(grid,cluster)
        if len(cluster)>0.99*N*N:
            count = count + 1
        if count>100:
            break
        if count ==1:
            index = i
        
        magnetizations.append(np.mean(grid))
    
    if mags:
        magnetizations = np.abs(magnetizations)
        return grid,np.mean(magnetizations[max(8*len(magnetizations)//9,index):])
    return grid

def N_behaviour_wolff(N_list,Klow,dK,steps,iters):
    """
    Obtains the behaviour of the magnetization as a function
    of the number of spins N in N_list
    """
    magnetizations = []
    magnetizations_std = []
    K_values = []
    
    for N in N_list:
        M = 1000
        grid = lattice(2,N)
        
        k_vals,mags,stds = ising2_magnetization(grid,Klow,dK,steps,N,M,iters)
        
        K_values.append(k_vals)
        magnetizations.append(mags)
        magnetizations_std.append(stds)
        
    return N_list, K_values, magnetizations,magnetizations_std
