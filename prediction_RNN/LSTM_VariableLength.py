# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 10:12:05 2017

@author: Administrator
"""


#%% import modules
import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt


#%% hyperparameters
BATCH_START = 0
MAX_STEPS = 20  ## number of residues, this should be variant in different batches
BATCH_SIZE = 50 ## number of proteins, this should be variant in different bathses
INPUT_SIZE = 2   ##
OUTPUT_SIZE = 2  ## angles
CELL_SIZE = 64
LR = 0.006
BATCH_START_TEST = 0


#%% ===========================================================================
# customize functions
#%%
# function to create test input and output sequence with variable length
def get_batch(timeSteps=MAX_STEPS, plotSwitch=False):
    global BATCH_START # declear as gloabal variables to modify their values inside function
    # xs shape (BATCH_SIZE, timeSteps)
    xs = np.arange(BATCH_START, BATCH_START+timeSteps*BATCH_SIZE).reshape((BATCH_SIZE, timeSteps)) / 10 * np.pi
    xs = xs[:,:,np.newaxis]

    # input (BATCH_SIZE, timeSteps, input_size)
    seq = np.append(np.sin(xs), np.cos(xs), axis=2)
    # output (BATCH_SIZE, timeSteps, output_size)
    res = np.append(np.sin(3*xs)**2+3*np.cos(xs+np.pi/6.0)**3, np.tanh(3*np.sin(xs)+np.cos(xs+np.pi/3)), axis=2)

    # increase the start of batch by timeSteps
    BATCH_START += timeSteps
    # plot res and seq in the first batch
    if plotSwitch:
        plt.subplot(211)
        plt.plot(xs[0,:], seq[0,:,0], 'r', xs[0,:], seq[0,:,1], 'b--')
        plt.ylabel('input')
        plt.subplot(212)
        plt.plot(xs[0,:], res[0,:,0], 'r', xs[0,:], res[0,:,1], 'b--')
        plt.ylabel('output')
        plt.show()
    # return seq, res and shape (size, step, input)
    return [seq, res, xs]


#%% length of sequence in the same batch
# assume that the sequences are padded with zero vectors to fill up the remaining time steps in the batch
def length(sequence):
    used = tf.sign(tf.reduce_max(tf.abs(sequence), reduction_indices=2)) # take sign of max abs along feature axis
    length = tf.reduce_sum(used, reduction_indices=1) # sum along time axis
    length = tf.cast(length, tf.int32) # data type -> int32
    return length


#%%
# define class for LSTMRNN for variable length sequence
class LSTMRNN(object):
    # initializer
    def __init__(self, max_steps, input_size, output_size, cell_size, batch_size):
        self.max_steps = max_steps
        self.input_size = input_size
        self.output_size = output_size
        self.cell_size = cell_size
        self.batch_size = batch_size
        # placeholders for inputs
        with tf.name_scope('inputs'):
            # placeholder for input: (batch_size, max_steps, input_size)
            self.xs = tf.placeholder(tf.float32, [None, max_steps, input_size], name='xs')
            # placeholder for output: (batch_size, max_steps, output_size)
            self.ys = tf.placeholder(tf.float32, [None, max_steps, output_size], name='ys')
        # variables for input hidden layer
        with tf.variable_scope('in_hidden'):
            self.add_input_layer()
        # variables for LSTM cell # NOTE: modified into multilayers!
        with tf.variable_scope('LSTM_cell'):
            self.add_cell()
        # variables for output hidden layer
        with tf.variable_scope('out_hidden'):
            self.add_output_layer()
        # cost function
        with tf.name_scope('cost'):
            self.compute_cost_length() #self.compute_cost()
        # train optimizer
        with tf.name_scope('train'):
            self.train_op = tf.train.AdamOptimizer(LR).minimize(self.cost)
    
    # add one input hidden layer
    def add_input_layer(self):
        # (batch_size, max_steps, in_size) ==> (batch_size*max_steps, in_size)
        l_in_x = tf.reshape(self.xs, [-1, self.input_size], name='2_2D')
        # Ws (in_size, cell_size)
        Ws_in = self._weight_variable([self.input_size, self.cell_size])
        # bs (cell_size,)
        bs_in = self._bias_variable([self.cell_size,])
        # l_in_y (batch_size*max_steps, cell_size)
        with tf.name_scope('Wx_plus_b'):
            l_in_y = tf.matmul(l_in_x, Ws_in) + bs_in
        # reshape l_in_y ==> (batch_size, max_steps, cell_size)
        self.l_in_y = tf.reshape(l_in_y, [-1, self.max_steps, self.cell_size], name='2_3D')
    
    # add one LSTM-RNN cell: (batch_size, max_steps, cell_size) -> (batch_size, max_steps, cell_size)
    def add_cell(self):
        # create basic LSTM cell --> http://arxiv.org/abs/1409.2329
        lstm_cell = tf.contrib.rnn.BasicLSTMCell(self.cell_size, forget_bias=1.0, state_is_tuple=True)
        # create initial state as tuple:(batch_size, cell_size, ?)
        with tf.name_scope('initial_state'):
            self.cell_init_state = lstm_cell.zero_state(self.batch_size, dtype=tf.float32) # NOTE:explicitly use batch_size
        # creates a recurrent neural network specified by RNNCell --> https://www.tensorflow.org/api_docs/python/tf/nn/dynamic_rnn
        self.cell_outputs, self.cell_final_state = tf.nn.dynamic_rnn(
                lstm_cell, self.l_in_y, initial_state=self.cell_init_state, time_major=False, sequence_length=length(self.xs))
    
    # add one output hidden layer
    def add_output_layer(self):
        # (batch_size, max_steps, cell_size) ==> (batch_size*max_steps, cell_size)
        l_out_x = tf.reshape(self.cell_outputs, [-1, self.cell_size], name='2_2D')
        # Ws (cell_size, out_size)
        Ws_out = self._weight_variable([self.cell_size, self.output_size])
        # bs (out_size,)
        bs_out = self._bias_variable([self.output_size,])
        # (batch*max_steps, out_size)
        with tf.name_scope('Wx_plus_b'):
            self.pred = tf.matmul(l_out_x, Ws_out) + bs_out
    
    # compute cost function by the actual length
    def compute_cost_length(self):
        # reshape prediction value to (batch_size, max_steps, out_size)
        pred_3D = tf.reshape(self.pred, [-1, self.max_steps, self.output_size], name='pred_3D')
        # compute losses for each step
        with tf.name_scope('losses'):
            # loss on each step (batch_size, max_steps)
            losses = tf.reduce_sum(self.ms_error(pred_3D, self.ys), reduction_indices=2)
            # mask of padding part
            mask = tf.sign(tf.reduce_max(tf.abs(self.ys), reduction_indices=2))
            # mean losses of each sequence: average along steps
            # shape: (batch_size,)
            losses = tf.div(tf.reduce_sum(losses * mask, reduction_indices=1), # total loss of each sequence by the actual length
                            tf.reduce_sum(mask, reduction_indices=1), # divide with actual length
                            name='average_losses_per_seq') 
        # average loses over batch
        with tf.name_scope('average_cost'):
            self.cost = tf.div(
                    tf.reduce_sum(losses, name='losses_sum'), # sum over batches
                    self.batch_size,
                    name='average_cost_per_batch')
        # record cost into summary
        tf.summary.scalar('cost', self.cost)
        
    # compute cost function
    def compute_cost(self):
        # loss on each step (batch_size*max_steps)
        losses = tf.contrib.legacy_seq2seq.sequence_loss_by_example(
                [tf.reshape(self.pred, [-1], name='reshape_pred')],
                [tf.reshape(self.ys, [-1], name='reshape_target')],
                [tf.ones([self.batch_size * self.max_steps * self.output_size], dtype=tf.float32)],
                average_across_timesteps=True,
                softmax_loss_function=self.ms_error,
                name='losses')
        # average loses over batch
        with tf.name_scope('average_cost'):
            self.cost = tf.div(
                    tf.reduce_sum(losses, name='losses_sum'),
                    self.batch_size,
                    name='average_cost')
        # record cost into summary
        tf.summary.scalar('cost', self.cost)
    
    # mean squared error
    def ms_error(self, y_pre, y_target):
        return tf.square(tf.subtract(y_target, y_pre))
    
    # weights: initialized with normal distribution
    def _weight_variable(self, shape, name='weights'):
        initializer = tf.random_normal_initializer(mean=0., stddev=1.,)
        # details of get_variable() --> https://www.tensorflow.org/programmers_guide/variable_scope
        return tf.get_variable(shape=shape, initializer=initializer, name=name)
    
    # biases: initialized with small positive constant
    def _bias_variable(self, shape, name='biases'):
        initializer = tf.constant_initializer(0.1)
        return tf.get_variable(shape=shape, initializer=initializer, name=name)


#%% ===========================================================================
#%%
if __name__ == '__main__':
    # create an instance of LSTMRNN
    model = LSTMRNN(MAX_STEPS, INPUT_SIZE, OUTPUT_SIZE, CELL_SIZE, BATCH_SIZE)
    
    # create a session
    sess = tf.Session()
    
    # for tensorboard 
    merged = tf.summary.merge_all()
    writer = tf.summary.FileWriter("logs_LSTM", sess.graph)
    # to see the graph in command line window, then type:
    #   python -m tensorflow.tensorboard --logdir=logs_Regression

    # initialze all variables
    init = tf.global_variables_initializer()
    sess.run(init)
    
    # open figure to plot
    plt.ion()
    plt.show()
    
    # total number of runs
    num_run = 200
    # number of time steps in each run
    steps = np.random.randint(MAX_STEPS//3, MAX_STEPS+1, num_run)
    
    for i in range(num_run):
        # obtain one batch
        seq, res, xs = get_batch(timeSteps=steps[i])
        # padding to max_steps
        seq_padding = np.append(seq, np.zeros([BATCH_SIZE, MAX_STEPS-steps[i], INPUT_SIZE]), axis=1)
        res_padding = np.append(res, np.zeros([BATCH_SIZE, MAX_STEPS-steps[i], OUTPUT_SIZE]), axis=1)
        
        # create the feed_dict
        if i == 0:
            feed_dict = {
                    model.xs:seq_padding,
                    model.ys:res_padding
                    } # the first input, no need to input initial_state
        else:
            feed_dict = {
                    model.xs:seq_padding,
                    model.ys:res_padding,
                    model.cell_init_state:state # update current state
                    }
        # run one step of training
        _, cost, state, pred = sess.run(
                [model.train_op, model.cost, model.cell_final_state, model.pred],
                feed_dict=feed_dict)
        # plotting
        plt.subplot(211)
        plt.plot(xs[0,:], res[0,:,0].flatten(), 'r', xs[0,:], pred[:,0].flatten()[:steps[i]], 'b--')
        plt.ylim((-4, 4))
        plt.ylabel('output_feature_1')
        plt.subplot(212)
        plt.plot(xs[0,:], res[0,:,1].flatten(), 'r', xs[0,:], pred[:,1].flatten()[:steps[i]], 'b--')
        plt.ylim((-2, 2))
        plt.ylabel('output_feature_2')
        plt.draw()
        plt.pause(0.3)
        # write to log
        if i % 20 == 0:
            print('cost: ', round(cost, 4))
            result = sess.run(merged, feed_dict)
            writer.add_summary(result, i)
        

