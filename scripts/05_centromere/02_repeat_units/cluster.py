import sys
from collections import defaultdict, deque

edges = defaultdict(set)
for line in open("self.sim90.tsv"):
    q,s,*_ = line.rstrip("\n").split("\t")
    edges[q].add(s)
    edges[s].add(q)

seen=set()
cid=0
for node in edges:
    if node in seen: 
        continue
    cid += 1
    dq=deque([node])
    seen.add(node)
    comp=[]
    while dq:
        x=dq.popleft()
        comp.append(x)
        for y in edges[x]:
            if y not in seen:
                seen.add(y)
                dq.append(y)
    print(f"Cluster{cid}\t" + "\t".join(sorted(comp)))
