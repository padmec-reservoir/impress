M.happiness[:] = -1
for x in range(len(M.coarse.interfaces_nodes)):
    M.happiness[M.coarse.interfaces_nodes[x]] = x



M.joy[:] = -1
for x in range(len(M.coarse.interfaces_edges)):
    M.joy[M.coarse.interfaces_edges[x]] = x


M.pride[:] = -1
for x in range(len(M.coarse.interfaces_faces)):
    M.pride[M.coarse.interfaces_faces[x]] = x


M.dreams[:] = -1
for x in range(len(M.coarse.interfaces_faces)):
    M.pride[M.coarse.interfaces_faces[x]] = x
