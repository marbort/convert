import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse
def rotate_necklace(necklace:list,stars:list,turns:int):
    for i in range(turns):
        
        necklace = [necklace * stars for necklace, stars in zip(necklace, stars)]
        necklace = [necklace[-1]] + necklace[:-1]
        return necklace

def create_necklace(num_beads:int, num_stars:int):
    rng=np.random.default_rng()
    black=rng.integers(low=0,high=num_beads,size=np.random.randint(1,num_beads))
    
    necklace = np.ones(num_beads)
    for i in black:
        necklace[i]=-1
        
    stars = np.ones(num_beads)
    try:
        active_stars = rng.integers(low=0,high=num_beads,size=np.random.randint(1,num_stars))
    except:
        active_stars = rng.integers(low=0,high=num_beads,size=1)
    for i in active_stars:
        stars[i]=-1
    
    assert num_stars < num_beads, "number of stars must be less than number of beads"
    
    return necklace, stars


parser=argparse.ArgumentParser(description="Simulate a necklace with stars")

parser.add_argument('--beads', type=int, help="number of beads in the necklace",default=10)
parser.add_argument('--stars', type=int, help="number of stars in the necklace",default=3)
parser.add_argument('--turns', type=int, help="number of turns to simulate",default=10)
parser.add_argument('--plot', action='store_true', help="plot the necklace")
args=parser.parse_args()

#necklace = [1,-1,-1,1,1,-1,-1]
#stars=[-1,1,1,1,-1,1,1]
colors=[]
turns=[]
necklace,stars=create_necklace(args.beads,args.stars)
colors.append(sum(necklace))
turns.append(0)

for i in range(1,args.turns):   
    necklace=rotate_necklace(necklace,stars,i)
    color=sum(necklace)
    colors.append(color)
    turns.append(i)
    #print(f"color after {i} turns: {color}")
print(np.mean(colors))
if args.plot:
    print(necklace,stars)
    plt.plot(turns,colors,marker='o')
    plt.title(f"Necklace with {args.beads} beads and {args.stars} Stars")
    plt.show()

    