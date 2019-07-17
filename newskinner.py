def skinner_operation(self):
    self.skin = sk.Skinner(self.mb)
    print("Entering skinner test")

    if self.dimension == 3:
        faces_on_skin_handles = self.skin.find_skin(self.root_set, self.all_volumes[:])
        edges_on_skin_handles = self.mtu.get_bridge_adjacencies(faces_on_skin_handles, 2, 1)
        nodes_on_skin_handles = self.mtu.get_bridge_adjacencies(faces_on_skin_handles, 2, 0)
        volumes_on_skin_handles = self.mtu.get_bridge_adjacencies(faces_on_skin_handles, 0, 3)
    elif self.dimension == 2:
        edges_on_skin_handles = skin.find_skin(self.root_set, self.all_faces[:])
        nodes_on_skin_handles = self.access_handle(edges_on_skin_handles)
        nodes_in_faces = ([self.mb.get_adjacencies(el_handle,0) for el_handle in self.all_faces])
        check_faces= np.asarray([rng.intersect(el_handle,nodes_on_skin_handles) for el_handle in nodes_in_faces])
        external_faces_index = np.array([el_handle.empty() for el_handle in check_faces]).astype(bool)
        faces_on_skin_handles = self.range_index(np.bitwise_not(external_faces_index),self.all_faces)
        volumes_on_skin_handles = rng.Range()

    print("Skinning Operation Successful")
    return [nodes_on_skin_handles, edges_on_skin_handles, faces_on_skin_handles, volumes_on_skin_handles]
