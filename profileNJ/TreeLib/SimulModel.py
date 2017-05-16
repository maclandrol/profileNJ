# Tree simulation
# part of this source code was inspirer from the dendropy project 
# the two concerned method are  the pure_birth and birth_death function
# with heavy modification on the birth_death function
#

from TreeClass import TreeClass
import random
from collections import defaultdict as ddict
import logging

INF = float('Inf')

def stop_with_tree_size(tree, **kwargs):
    nsize = kwargs['nsize']
    cur_size = kwargs['cur_size']
    return nsize == cur_size

def stop_with_max_time(tree, **kwargs):
    ctime = kwargs['cur_time']
    mtime = kwargs['max_time']
    return mtime <= ctime

def stop_with_size_range(tree, **kwargs):
    size_range = kwargs.get('size_range', None)
    cur_size = kwargs['cur_size']
    try:
        if len(size_range) == 2:
            return size_range[0] <= cur_size and cur_size <= size_range[1]
    except:
        return False
    


class TotalExtinction(Exception):
    """Exception to be raised when branching process results in all lineages going extinct."""
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)


class MissingParameterError(Exception):
    """Exception to be raised when parameters are missing."""
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)



class FunctionSlot(object):
    """Light functionSlot class
    """

    def __init__(self, name="Anonymous Function", rand_apply=False):
        """ The creator of the FunctionSlot Class """
        self.funcList = []
        self.slotName = name

    def __typeCheck(self, func):
        """ Used internally to check if a function passed to the
        function slot is callable. Otherwise raises a TypeError exception.

        :param func: the function object
        """
        if not callable(func):
            Util.raiseException("The function must be a method or function", TypeError)

    def __iadd__(self, func):
        """ To add more functions using the += operator

            .. versionadded:: 0.6
                The __iadd__ method.
        """
        self.__typeCheck(func)
        self.funcList.append(func)
        return self

    def __getitem__(self, index):
        """ Used to retrieve some slot function index """
        return self.funcList[index]

    def __setitem__(self, index, value):
        """ Used to set the index slot function """
        self.__typeCheck(value)
        self.funcList[index] = value

    def __iter__(self):
        """ Return the function list iterator """
        return iter(self.funcList)

    def __len__(self):
        """ Return the number of functions on the slot

        .. versionadded:: 0.6
            The *__len__* method
        """
        return len(self.funcList)

    def clear(self):
        """ Used to clear the functions in the slot """
        if len(self.funcList) > 0:
            del self.funcList[:]

    def add(self, func):
        """ Used to add a function to the slot

        :param func: the function to be added in the slot
        """
        self.__typeCheck(func)
        self.funcList.append(func)

    def isEmpty(self):
        """ Return true if the function slot is empy """
        return (len(self.funcList) == 0)

    def set(self, func):
        """ Used to clear all functions in the slot and add one

        :param func: the function to be added in the slot
        .. note:: the method *set* of the function slot remove all previous
                     functions added to the slot.
        """
        self.clear()
        self.__typeCheck(func)
        self.add(func)

    def apply(self, index, obj, **args):
        """ Apply the index function

        :param index: the index of the function
        :param obj: this object is passes as parameter to the function
        :param args: this args dictionary is passed to the function

        """
        if len(self.funcList) <= 0:
            raise Exception("No function defined: for " + self.slotName)
        return self.funcList[index](obj, **args)

    def applyFunctions(self, obj=None, **args):
        """ Generator to apply all function slots in obj

        :param obj: this object is passes as parameter to the function
        :param args: this args dictionary is passed to the function

        """
        if len(self.funcList) <= 0:
            Util.raiseException("No function defined: " + self.slotName)

        for f in self.funcList:
            yield f(obj, **args)

    def __repr__(self):
        """ String representation of FunctionSlot """
        strRet = "Slot [%s] (Count: %d)\n" % (self.slotName, len(self.funcList))

        if len(self.funcList) <= 0:
            strRet += "\t\tNo function\n"
            return strRet

        for f in self.funcList:
            fname = "-"
            try: 
                fname = f.func_name
            except AttributeError:
                fname =  f.func.func_name
            strRet += "\t\tName: %s\n" % (fname)

        return strRet


class SimulModel(object):
    def __init__(self, stopcrit=None, seed=None, debug=False):
        self.debug = debug
        if stopcrit is None:
            stopcrit = FunctionSlot("Stopping Crit")
        self.stopcrit = stopcrit
        random.seed(seed)   
        logging.basicConfig(level=logging.DEBUG)#, format='%(relativeCreated)6d %(threadName)s %(message)s')


    def add_stopping_crit(self, func):
        self.stopcrit.add(func)

    def set_stopping_crit(self, func):
        self.stopcrit.set(func)

    def clear_stopping_crit(self):
        self.stopcrit.clear()

    def debug_msg(self, msg):
        if self.debug:
            logging.debug(msg)


    def pure_birth_tree(self, birth=1.0, **kwargs):
        """Generates a uniform-rate pure-birth process tree.
        You can pass  supplemental argument:
            - ``nsize`` : total number of leaves before stopping
            - ``names_library`` : list of names for the leaves
            - ``max_time`` : maximum time for simulation
            - if nsize if not provided, the length of names_library will be used
            - if nsize is larger than ``names_library``, leaves name will be completed with 
            random names in the following format : "T%d" (T1, T2, etc)
        and 
        """
        tree = TreeClass()
        tree.dist = 0.0
        # time of waiting
        # compared to parent
        done = False
        tname = kwargs.get("names_library", [])
        nsize = kwargs.get("nsize", len(tname))
        max_time = kwargs.get("max_time", None)
        pb_stop = FunctionSlot("Pure birth stopping")
        if nsize:
            pb_stop.add(stop_with_tree_size)
        if max_time:
            pb_stop.add(stop_with_max_time)
        if pb_stop.isEmpty() and self.stopcrit.isEmpty():
            raise MissingParameterError("Either specify a names_library, nsize, max_time or a stopping criterion")
        
        extra_param = {}
        for k,v in kwargs.items():
            if k not in ['nsize', 'max_time', 'removeloss']:
                extra_param[k] = v

        extra_param['nsize'] =  nsize
        extra_param['max_time'] =  max_time


        # fill namespace to desired size
        total_time = 0
        while True:
            # time before new node
            # given the probability of birth
            leaf_nodes = tree.get_leaves()
            wtime = random.expovariate(len(leaf_nodes)/birth)
            total_time += wtime
            for leaf in leaf_nodes:
                leaf.dist += wtime
            if not pb_stop.isEmpty():
                for val in pb_stop.applyFunctions(tree, cur_time=total_time, cur_size=len(leaf_nodes), **extra_param):
                    done = done or val
            if not self.stopcrit.isEmpty():
                for val in self.stopcrit.applyFunctions(tree,  cur_time=total_time, cur_size=len(leaf_nodes), **extra_param):
                    done = done or val
            if done:
                break

            if max_time is None or total_time <= max_time:
                # now add new node to a random leaf
                node = random.choice(leaf_nodes)
                c1 = TreeClass()
                c2 = TreeClass()
                node.add_child(c1)
                node.add_child(c2)
                c1.dist = 0.0
                c2.dist = 0.0

        leaf_nodes = tree.get_leaves()
        leaf_compteur = 1
        total_time += wtime
        for ind, node in enumerate(leaf_nodes):
            if ind < len(tname):
                node.name = tname[ind]
            else:
                node.name = "T%d"%leaf_compteur
                leaf_compteur += 1
        return tree


    def birth_death_tree(self, birth, death, **kwargs):
        """
        Returns a birth-death tree with birth rate specified by ``birth``, and
        death rate specified by ``death``, and  edge lengths in continuous (real)
        units.

        You can pass  supplemental argument:
        - ``nsize`` : total number of leaves before stopping
        - ``names_library`` : list of names for the leaves
        - ``max_time`` : maximum time for simulation
        - if nsize if not provided, the length of names_library will be used
        - if nsize is larger than ``names_library``, leaves name will be completed with 
        random names in the following format : "T%d" (T1, T2, etc)
        - If `max_time` is given as a keyword argument, tree is grown for
        a maximum of ``max_time``.
        - if `removeloss` is given as argument (default True), extinct taxa are removed

        Under some conditions, it is possible for all lineages on a tree to go extinct.
        In this case, if the keyword argument ``repeat_until_success`` is |True| (default),
        then a new branching process is initiated. Otherwise a TotalExtinction error is raised.

        """

        tree = TreeClass()
        tree.dist = 0.0

        done = False
        removeloss = kwargs.get("removeloss", True)
        repeat_until_success = kwargs.get("repeat_until_success", True)
        names_library = kwargs.get("names_library", [])
        nsize = kwargs.get("nsize", len(names_library))
        max_time = kwargs.get("max_time", None)
        pb_stop = FunctionSlot("birth death stopping")
        if nsize:
            pb_stop.add(stop_with_tree_size)
        if max_time:
            pb_stop.add(stop_with_max_time)
        if pb_stop.isEmpty() and self.stopcrit.isEmpty():
            raise MissingParameterError("Either specify a names_library, nsize, max_time or a stopping criterion")
        
        extra_param = {}
        for k,v in kwargs.items():
            if k not in ['nsize', 'max_time', 'removeloss']:
                extra_param[k] = v

        extra_param['nsize'] =  nsize
        extra_param['max_time'] =  max_time
        
        # initialize tree
        tree = TreeClass()
        tree.dist = 0.0

        #_LOG.debug("Will generate a tree with no more than %s leaves to get a tree of %s leaves" % (str(gsa_ntax), str(nsize)))
        leaf_nodes = tree.get_leaves()
        curr_num_leaves = len(leaf_nodes)
          
        total_time = 0

        died = set([])
        event_rate = float(birth + death)  

        while True:
            # waiting time based on event_rate
            wtime = random.expovariate(event_rate)
            #_LOG.debug("Drew waiting time of %f from hazard parameter of %f" % (wtime, all_rates))

            total_time += wtime
            for leaf in leaf_nodes:
                # extinct leaves cannot update their branches length
                if not leaf.has_feature('name', name=TreeClass.LOST):
                    leaf.dist += wtime

            if not pb_stop.isEmpty():
                for val in pb_stop.applyFunctions(tree, cur_time=total_time, cur_size=curr_num_leaves, **extra_param):
                    done = done or val
            if not self.stopcrit.isEmpty():
                for val in self.stopcrit.applyFunctions(tree,  cur_time=total_time, cur_size=curr_num_leaves, **extra_param):
                    done = done or val
            if done:
                break
            # if event occurs within time constraints
            if max_time is None or total_time <= max_time:

                # select node at random, then find chance it died or give birth (speciation)
                node = random.choice(leaf_nodes)
                eprob = random.random()
                leaf_nodes.remove(node)
                curr_num_leaves -= 1
                if eprob < birth/event_rate:
                    #_LOG.debug("Speciation")
                    c1 = TreeClass()
                    c2 = TreeClass()
                    c1.dist = 0
                    c2.dist = 0
                    node.add_features(type=TreeClass.SPEC)
                    node.add_child(c1)
                    node.add_child(c2)
                    leaf_nodes.append(c1)
                    leaf_nodes.append(c2)
                    curr_num_leaves += 2
                else:
                    #_LOG.debug("Extinction")
                    if curr_num_leaves > 0:
                        #_LOG.debug("Will delete " + str(id(nd)) + " with parent = " + str(id(nd.parent_node)))
                        died.add(node)
                        node.add_features(type=TreeClass.LOST)
                    else:
                        if not repeat_until_success:
                            raise TotalExtinction("All lineage went extinct, please retry")
                        # Restart the simulation because the tree has gone extinct
                        tree = TreeClass()
                        leaf_nodes = tree.get_leaves()
                        curr_num_leaves = 1
                        died = set([])
                        total_time = 0
                
                # this will always hold true
                assert curr_num_leaves == len(leaf_nodes)

        if removeloss:
            leaves = set(tree.get_leaves()) - died
            tree.prune(leaves)
            tree.delete_single_child_internal(enable_root=True)

        leaf_nodes = tree.get_leaves()
        #wtime = random.expovariate(event_rate)
        leaf_compteur = 1
        nlc = 0
        for ind, node in enumerate(leaf_nodes):
            if not node.has_feature('type', name=TreeClass.LOST):
                #node.dist += wtime
                if nlc < len(names_library):
                    node.name = names_library[nlc]
                    nlc += 1
                else:
                    node.name = "T%d"%leaf_compteur
                    leaf_compteur += 1
        return tree


    def dual_gtree_sptree(self):
        pass


    def dlt_tree_from_sptree(self, sptree, birth, death, transfer=0.0, **kwargs):
        """Simulate a gene tree within a species tree with birth death and transfer rates
            - ``removeloss`` set to True to remove losses from the simulated tree (default: True)
            - ``disallow_suc_trn`` set to True to disable successive transfer. If a node parent has a transfer event,
            the node can't have a transfer event too. This prevent circular event
            - ``names_library`` : list of gene label at the leaves. If there is not enough names, they will be 
            duplicated (ex : for 5 leaves and [a, b], leaves name will be [a, b, aa, bb, aaa]
            - `nsize``: number of non-extinct leaves wished (almost useless if repeat_until_success is false)
            - ``repeat_until_success`` : repeat simulation until `nsize` leaves are obtained or stopping criterion is attained
            if ``repeat_until_success`` is True and there are not stopping criterion (nsize of model.stopcrit), simulation will run until
            the first tree where all lineage are not extinct is obtaine


        """
        repeat_until_success = kwargs.get("repeat_until_success", True)
        # sanity check
        if len(sptree) < 2:
            raise ValueError("Species Tree should have at least 2 leaves")
        nsize = kwargs.get('nsize', None)
        ntry = 0
        done = False
        last_was_extinct = False
        no_stopping = self.stopcrit.isEmpty() and not nsize
        while True:
            extra_param = {}
            try:
                gtree, recon, events, ecounter, transfers = self.sample_from_tree(sptree, birth, death, transfer, **kwargs)
                cur_size = len(gtree.get_leaves(is_leaf_fn=lambda x: x.is_leaf() and not x.has_feature('type', TreeClass.LOST)))
                if nsize and cur_size == nsize:
                    done = True
                # all leaves are extincts
                if cur_size == 0:
                    raise TotalExtinction
                if not self.stopcrit.isEmpty():
                    for val in self.stopcrit.applyFunctions(gtree, cur_size=cur_size, **kwargs):
                        done = done or val

                last_was_extinct = False
            except TotalExtinction:
                last_was_extinct = True

            # successfully stopped
            # or no repeat requested
            # or no stopping crit and not extinct
            if repeat_until_success and no_stopping and not last_was_extinct:
                break
            elif not repeat_until_success or done:
                break
            ntry += 1

        if last_was_extinct:
            raise TotalExtinction("All lineage went extinct during simulation")
        
        event_logger = {}
        event_logger['recon'] = recon
        event_logger['count'] = ecounter
        event_logger['transfers'] = transfers
        event_logger['events'] = events

        return gtree, event_logger


    def sample_from_tree(self, sptree, birth, death, gain, **kwargs):
        """Sample a tree within another tree using the rate specified 
            Note that a tree with all leaves being extinct can be returned by this function
            Use dlt_tree_from_sptree if you want to prevent this.
        """
        
        # initialize gene tree
        sptree.compute_branches_length()
        sptree.label_internal_node()
        removeloss = kwargs.get("removeloss", True)
        disallow_suc_trn = kwargs.get("disallow_suc_trn", True)
        leave_names = kwargs.get("names_library", [])

        gtree = TreeClass()
        recon = {gtree: sptree}
        events = {gtree: "spec"}
        losses = set()
        transfers = {}
        true_event_counter=ddict(int)
        snode_counter = ddict(int)
        if not leave_names:
            leave_names = lambda sp, x: sp + "_" + str(x)
        name_counter = 0
        def create_history(snode, gnode):
            if snode.is_leaf():
                if isinstance(leave_names, list):
                    n_encounter = name_counter/len(leave_names)
                    gnode.name = leave_names[name_counter%len(leaves_name)]
                    gnode.name += ("_"+gnode.name)*n_encounter
                else:
                    snode_counter[snode.name] += 1
                    gnode.name = leave_names(snode.name, snode_counter[snode.name])
                events[gnode] = "leaf"
                gnode.add_features(type=TreeClass.SPEC)
            else:
                for schild in snode.get_children():
                    # get branches event for branch (snode, schild)
                    # during time = schild.dist
                    recnode, died, transfered, smap = self.sample_event_on_branches(schild.dist,
                        schild, birth, death, gain, keeplosses=(not removeloss), ign_suc_trn=disallow_suc_trn, ecounter=true_event_counter)
                    gnode.add_child(recnode)
                    # update ist of losses
                    losses.update(died)
                    transfers.update(transfered)
                    recon.update(smap)
                    next_cand = []
                    # then record reconciliation that happened
                    #print recnode.get_ascii(attributes=[], show_internal=True)
                    #print schild
                    for node in recnode.traverse():
                        node.add_features(species=recon[node].name)
                        if node.type == TreeClass.LOST:
                            events[node] = "loss"
                            # died at the start of coalescence
                            losses.add(node)
                        elif node.is_leaf():
                            node.add_features(type=TreeClass.SPEC)
                            events[node] = "spec"
                            next_cand.append(node)
                        elif node.type == TreeClass.AD:
                            events[node] = "dup"
                        else:
                            events[node] = "transfer"
                    
                    for new_node in next_cand:
                        create_history(recon[new_node], new_node)
                
                # if no child for node then it is a loss
                if gnode.is_leaf():

                    losses.add(gnode)
        create_history(sptree, gtree)
        
        gtree.delete_single_child_internal(enable_root=True)
        if removeloss:
            for node in gtree.traverse():
                if node in losses:
                    node.delete()
            gtree.delete_single_child_internal(enable_root=True)
            remove_from_history = set(recon.keys()) - set(gtree.traverse())
            for node in remove_from_history:
                del recon[node]
                if node in events.keys():
                    del events[node]

        if len(gtree) <= 1:
            raise TotalExtinction("All taxa are extinct.")        
        return gtree, recon, events, true_event_counter, transfers


    def sample_event_on_branches(self, time, spnode, birth, death, transfer, gnode=None, keeplosses=False, ign_suc_trn=False, ecounter={}):
        """Simulate a reconstructed birth death tree"""
        
        # we are going with a poisson process 
        # so the rate of having an event is
        # just the sum of rate
        event_rate = float(birth + death + transfer)
        died = set()
        transfered = {}
        map_to_spec = {}
        # create starting node if one is not given
        if gnode is None:
            gnode = TreeClass()
            map_to_spec[gnode] = spnode

        def event_in_time(time, node, spnode):
            # time for an event
            if event_rate == 0.0:
                next_t = INF
            else:
                next_t = random.expovariate(event_rate)
            
            if next_t > time:
                # no event on branch
                # we can stop
                node.dist = time
                node.add_features(type=INF)
            
            else:
                eprob = random.random()
                node.dist = next_t
                if eprob < birth*1.0 / event_rate:
                    # birth ==> duplication event
                    cnode = TreeClass()
                    node.add_child(cnode)
                    map_to_spec[cnode] = spnode
                    event_in_time(time - next_t, cnode, spnode)
                    # compute event on the remaining time
                    cnode = TreeClass()
                    node.add_child(cnode)
                    map_to_spec[cnode] = spnode
                    event_in_time(time - next_t, cnode, spnode)
                    node.add_features(type=TreeClass.AD)
                    ecounter['dup'] += 1
                
                elif eprob < (birth+death)*1.0/event_rate:
                    # death happen ==> loss
                    node.add_features(type=TreeClass.LOST)
                    map_to_spec[node] = spnode
                    ecounter['loss'] += 1
                    died.add(node)
                else:
                    # give gene to another species ==> transfer
                    contemp_transfer_nodes =  list(spnode.get_incomparable_list(timeconsistent=True, wtime=next_t))
                    if contemp_transfer_nodes and not(ign_suc_trn and node.up and node.up.has_feature('type', name=TreeClass.TRANSFER)):
                        cand_receiver = random.choice(contemp_transfer_nodes)
                       
                        node.add_features(type=TreeClass.TRANSFER)
                        ecounter['transfer'] += 1
                        cnode = TreeClass()
                        node.add_child(cnode)
                        map_to_spec[cnode] = spnode
                        event_in_time(time - next_t, cnode, spnode)
                        
                        cnode = TreeClass()
                        node.add_child(cnode)
                        cnode.add_features(transfered=True)
                        transfered[cnode] = cand_receiver
                        map_to_spec[cnode] = cand_receiver
                        t = cand_receiver.brlen - spnode.brlen + time - next_t
                        event_in_time(t, cnode, cand_receiver)
                    else:
                        # keep node as it is
                        # so speciation
                        node.add_features(type=INF)

                        self.debug_msg("Could not perform transfer at node")
                        self.debug_msg(node)

        event_in_time(time, gnode, spnode)

        if not keeplosses:
            leaves = set(gnode.get_leaves()) - died
      
            if len(leaves) == 0:
                gnode.add_features(type=TreeClass.LOST)
                died.add(gnode)
            else:
                gnode.prune(leaves)
            gnode.delete_single_child_internal()

        
        return gnode, died, transfered, map_to_spec


if __name__ == '__main__':
    # this are for test
    model = SimulModel(stopcrit=None, debug=True)
    #def pure_birth_tree(self, nsize, birth=1.0, leave_names=[]):

    sptree = model.pure_birth_tree(birth=0.5, nsize=5, names_library="abcdef")
    model.add_stopping_crit(stop_with_size_range)
    #sptree = model.birth_death_tree(birth=0.5, death=0.2, nsize=6, removeloss=False, names_library="abcdef", repeat_until_success=False)
    l = model.dlt_tree_from_sptree(sptree, 0.3, 0.4, removeloss=True, size_range=(7, 10), disallow_suc_trn=True, repeat_until_success=True)
    sptree.compute_branches_length()

    print sptree.get_ascii(show_internal=True, attributes=['name', 'dist'])
    gtree = l[0]
    print gtree.get_ascii(show_internal=True, attributes=['type', 'name', 'species'])