##  CPPN functions and classes

#   Imports
import math
import numpy

from evolve_soft_2d import utility

################################################################################

class cppn:
    """The CPPN class object
    """

    def __init__(
        self,
        seed: int,
        mod_n: int,
        scale: float,
        hl_n: int,
        hl_s: int,
        thresh: float,
        x: int,
        y: int,
        ) -> None:
        """The CPPN parameters

        Parameters
        ----------
        seed : int
            The seed for the random generation
        mod_n : int
            The number of models to be generated from a particular seed
        scale : float
            The scale of the focus on the model
        hl_n : int
            The number of hidden layers
        hl_s : int
            The size of the initial hidden layer
        thresh : float
            The rounding/removal threshold
        x : int
            The number of elements in the x-direction
        y : int
            The number of elements in the y-direction
        """

        self.seed = seed
        self.mod_n = mod_n
        self.scale = scale
        self.hl_n = hl_n
        self.hl_s = hl_s
        self.thresh = thresh
        self.x = x
        self.y = y
        
        #   The resolution of the grid
        self.res = self.x*self.y

        #   Build the grid
        self.grid = self.cppn_grid()

    def __repr__(self) -> str:
        """Format a representation of the CPPN

        Returns
        -------
        str
            Formatted representation of the CPPN for the log
        """

        r = "Model Dimensions:               {}x{} elements\n".format(self.x, self.y)
        r += "Model Seed:                     {}\n".format(self.seed)
        r += "Number Of Models Generated:     {}\n".format(self.mod_n)
        r += "Model Scale:                    1:{}\n".format(self.scale)
        r += "Number Of Hidden Layers:        {}\n".format(self.hl_n)
        r += "Size Of Initial Hidden Layer:   {}\n".format(self.hl_s)
        if self.thresh < 1:
            r += "Rounding Threshold:             {}\n".format(self.thresh)
        else:
            r += "Percentage Of Elements Removed: {}%\n".format(self.thresh)
        r += "Activation Functions:\n"
        for i in self.af:
            r += "{}\n".format(i)

        return r

    def cppn_grid(self) -> numpy.array:
        """Generates model grids

        Returns
        -------
        numpy.array
            The model grid
        """

        #   Initialisations
        self.af = []

        #   The list of possible activation functions
        af_l = [self.cppn_sin, self.cppn_cos, self.cppn_tanh, self.cppn_sigm, self.cppn_srel]
        af_o = [self.cppn_sigm, self.cppn_srel]

        #   Set the random generation seed
        numpy.random.seed(seed = self.seed)

        #   Generate the initial hidden layer for each model
        hl = numpy.random.uniform(low = -1, high = 1, size = (self.mod_n, self.hl_s)).astype(numpy.float32)

        #   Generate the grid matrix
        x_r = numpy.linspace(-1*self.scale, self.scale, num = self.x)
        x_m = numpy.matmul(numpy.ones((self.y, 1)), x_r.reshape((1, self.x)))

        y_r = numpy.linspace(-1*self.scale, self.scale, num = self.y)
        y_m = numpy.matmul(y_r.reshape((self.y, 1)), numpy.ones((1, self.x)))

        r_m = numpy.sqrt(x_m*x_m + y_m*y_m)

        x_d = numpy.tile(x_m.flatten(), self.mod_n).reshape(self.mod_n, self.res, 1)
        y_d = numpy.tile(y_m.flatten(), self.mod_n).reshape(self.mod_n, self.res, 1)
        r_d = numpy.tile(r_m.flatten(), self.mod_n).reshape(self.mod_n, self.res, 1)

        #   Scale the initial hidden layers
        hl_scale = numpy.reshape(hl, (self.mod_n, 1, self.hl_s))*numpy.ones((self.res, 1), dtype = numpy.float32)*self.scale

        #   Unwrap the grid matrices
        x_d_unwrap = numpy.reshape(x_d, (self.mod_n*self.res, 1))
        y_d_unwrap = numpy.reshape(y_d, (self.mod_n*self.res, 1))
        r_d_unwrap = numpy.reshape(r_d, (self.mod_n*self.res, 1))
        hl_unwrap = numpy.reshape(hl_scale, (self.mod_n*self.res, self.hl_s))

        #   Build the network
        n = self.fully_connected(hl_unwrap, self.hl_n, True, self.seed) + self.fully_connected(x_d_unwrap, self.hl_n, False, self.seed + 1) + self.fully_connected(y_d_unwrap, self.hl_n, False, self.seed + 2) + self.fully_connected(r_d_unwrap, self.hl_n, False, self.seed + 3)

        #   Transpose the network
        n = n.T

        if self.hl_n > 1:

            #   Loop through the second to second-last hidden layers
            for i in range(1, self.hl_n - 1):

                #   Set the seed for each layer
                numpy.random.seed(seed = self.seed + i)

                #   Select and record the activation function
                n[i], af_c = numpy.random.choice(af_l)(n[i - 1])
                self.af.append(af_c)

            #   Set the seed for the final layer
            numpy.random.seed(seed = self.seed)

            #   Apply and record the final function
            n[-1], af_o = numpy.random.choice(af_o)(n[-2])
            self.af.append(af_o)

        else:

            #   Set the seed for each layer
            numpy.random.seed(seed = self.seed)

            #   Select and record the activation function
            n[0], af_c = numpy.random.choice(af_l)(n[0])
            self.af.append(af_c)

            #   Apply and record the final function
            n[0], af_o = numpy.random.choice(af_o)(n[0])
            self.af.append(af_o)

        #   Reshape the grid to fit the given dimensions
        mod = numpy.reshape(n[-1], (self.mod_n, self.x, self.y))

        return mod

    def fully_connected(
        self,
        i_v: numpy.array,
        o_d,
        w_bias: bool,
        seed: int,
        ) -> numpy.array:
        """Connect all layers of the CPPN

        Parameters
        ----------
        i_v : numpy.array
            The input vector
        o_d
            The output dimensions
        seed : int
            The random generation
        w_bias : bool
            If the layers should be connected with bias

        Returns
        -------
        numpy.array
            The connected results
        """

        #   Set the random generation seed
        numpy.random.seed(seed = seed)

        #   Generate the random matrix
        m = numpy.random.standard_normal(size = (i_v.shape[1], o_d)).astype(numpy.float32)

        #   Multiply the input with the matrix
        result = numpy.matmul(i_v, m)

        #   Check if the bias must be included
        if w_bias:

            #   Generate the random bias
            bias = numpy.random.standard_normal(size = (1, o_d)).astype(numpy.float32)

            #   Add the bias to the result
            result += bias*numpy.ones((i_v.shape[0], 1), dtype = numpy.float32)

        return result

    def partly_connected(
        self,
        i_v: numpy.array,
        o_d,
        w_bias: bool,
        seed: int,
        ) -> numpy.array:
        """Connect a single layer of the hidden network

        Parameters
        ----------
        i_v : numpy.array
            The input vector
        o_d
            The dimensions of the output
        w_bias : bool
            The random generation
        seed : int
            If the layers should be connected with bias

        Returns
        -------
        numpy.array
            The connected results
        """

        #   Set the random generation seed
        numpy.random.seed(seed = seed)

        #   Generate the random matrix
        m = numpy.random.standard_normal(size = (i_v.shape[0], o_d)).astype(numpy.float32)

        #   Multiply the input with the matrix
        result = numpy.matmul(i_v, m)

        #   Check if the bias must be included
        if w_bias:

            #   Generate the random bias
            bias = numpy.random.standard_normal(size = (o_d)).astype(numpy.float32)
            
            #   Add the bias to the result
            result += bias.T

        return result

    def cppn_sin(
        self,
        hl: numpy.array
        ) -> (numpy.array, str):    
        """Apply sin as the activation function for the current layer

        Parameters
        ----------
        hl : numpy.array
            The current layer

        Returns
        -------
        numpy.array, str:
            The new layer
            The label of the activation function
        """

        name = "sin"

        out = numpy.sin(self.partly_connected(hl, self.res*self.mod_n, True, self.seed))

        return out, name

    def cppn_cos(
        self,
        hl: numpy.array
        ) -> (numpy.array, str):
        """Apply cos as the activation function for the current layer

        Parameters
        ----------
        hl : numpy.array
            The current layer

        Returns
        -------
        numpy.array, str:
            The new layer
            The label of the activation function
        """

        name = "cos"

        out = numpy.cos(self.partly_connected(hl, self.res*self.mod_n, True, self.seed))

        return out, name

    def cppn_tanh(
        self,
        hl: numpy.array
        ) -> (numpy.array, str):
        """Apply tanh as the activation function for the current layer

        Parameters
        ----------
        hl : numpy.array
            The current layer

        Returns
        -------
        numpy.array, str:
            The new layer
            The label of the activation function
        """

        name = "tanh"

        out = numpy.tanh(self.partly_connected(hl, self.res*self.mod_n, True, self.seed))

        return out, name

    def cppn_sigm(
        self,
        hl: numpy.array
        ) -> (numpy.array, str):
        """Apply a sigmoid as the activation function for the current layer

        Parameters
        ----------
        hl : numpy.array
            The current layer

        Returns
        -------
        numpy.array, str:
            The new layer
            The label of the activation function
        """

        name = "sigmoid"

        out = utility.sigmoid(self.partly_connected(hl, self.res*self.mod_n, True, self.seed))

        return out, name

    def cppn_srel(
        self,
        hl: numpy.array
        ) -> (numpy.array, str):
        """Apply smooth ReLu as the activation function for the current layer

        Parameters
        ----------
        hl : numpy.array
            The current layer

        Returns
        -------
        numpy.array, str:
            The new layer
            The label of the activation function
        """

        name = "smooth ReLu"

        out = utility.smooth_relu(self.partly_connected(hl, self.res*self.mod_n, True, self.seed))

        return out, name

################################################################################

class cppn_i:
    """The CPPN model
    """    

    def __init__(
        self,
        cppn: cppn,
        mod_id: int,
        ) -> None:
        """The CPPN model parameters

        Parameters
        ----------
        cppn : cppn
            The CPPN
        mod_id : int
            The model number
        """        

        self.cppn = cppn
        self.mod_id = mod_id

        #   The model grid
        self.grid = self.rem_thresh(self.cppn.grid[self.mod_id])

    def __repr__(self) -> str:
        """Format a representation of the CPPN model

        Returns
        -------
        str
            Formatted representation of the CPPN model for the log
        """

        r = "Model ID: {}\n".format(self.mod_id)
        r += "CPPN Parameters:\n{}".format(self.cppn)

        return r

    def rem_thresh(
        self,
        grid: numpy.array
        ) -> numpy.array:
        """Removes elements from a grid according to the specified threshold

        Parameters
        ----------
        grid : numpy.array
            The grid from which elements are to be removed

        Returns
        -------
        numpy.array
            The grid with the elements removed
        """

        #   Check if the threshold indicates that a percentage of elements should be removed
        if self.cppn.thresh > 1:

            #   Calculate the percentage as a decimal
            perc = self.cppn.thresh/100

            #   Calculate the number of elements to be removed
            b = int(math.ceil(self.cppn.x*self.cppn.y*perc))

            if b == self.cppn.x*self.cppn.y:

                #   Reshape the grid to be one-dimensional
                grid = numpy.zeros((self.cppn.x, self.cppn.y))

            else:

                #   Obtain the IDs of the elements to be removed
                ids = numpy.argpartition(grid, b, axis = None)

                #   Reshape the grid to be one-dimensional
                grid = grid.flatten()

                #   Loop through all elements
                for i in range(0, self.cppn.res):

                    #   Check if the current element ID refers to an element that needs to be removed
                    if i in ids[:b]:

                        #   Remove the element
                        grid[i] = 0

                    else:

                        #   Add the element
                        grid[i] = 1

                #   Reshape the grid to its original dimensions
                grid = grid.reshape(self.cppn.x, self.cppn.y)

        else:

            #   Loop through the grid
            for i in grid:

                for j in range(0, len(i)):

                    #   Check if the current value is above the threshold
                    if i[j] > self.cppn.thresh:

                        #   Set the element to exist
                        i[j] = 1

                    else:

                        #   Remove the element
                        i[j] = 0

        return grid

################################################################################

def cppn_rem(
    template,
    grid: list,
    ) -> list:
    """Obtain the list of elements to be removed according to the generated CPPN

    Parameters
    ----------
    template
        The unit template parameters
    grid : list
        The grid generated by the CPPN

    Returns
    -------
    list
        The list of element IDs to be removed
    """

    #   Initialisations
    rem = []

    #   Create a copy of the grid
    g_temp = grid[:]

    #   Reverse the order of the grid
    g_temp = g_temp[::-1]

    #   Calculate the element ID offset according to the boundary thickness
    offset = template.x_e*template.b + template.b + 1

    #   Loop through rows of the grid
    for i in g_temp:

        #   Loop through the elements in the current row
        for j in range(0, len(i)):

            #   Check if the current element needs to be removed
            if i[j] == 0:

                #   Add the element ID to the list of element IDs to be removed
                rem.append(j + offset)

        #   Increment the offset
        offset += template.x_e

    return rem