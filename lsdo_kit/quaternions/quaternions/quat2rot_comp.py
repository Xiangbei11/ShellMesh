import numpy as np

from openmdao.api import ExplicitComponent


class Quat2RotComp(ExplicitComponent):
    def initialize(self):
        self.options.declare('shape', types=tuple)
        self.options.declare('val', types=np.ndarray)
        self.options.declare('in_name',
                             types=str,
                             desc="This the name of the input quaternion")
        self.options.declare('out_name',
                             types=str,
                             desc="The name of the 3x3 rotation matrix")

    def setup(self):
        shape = self.options['shape']
        val = self.options['val']
        in_name = self.options['in_name']
        out_name = self.options['out_name']

        output_shape = (3, 3)  # The quaternion produces a 3x3 rotation matrix

        self.add_input(in_name, shape=shape, val=val)
        self.add_output(out_name, shape=output_shape)

        rows = np.tile(np.arange(np.prod(output_shape)), 4)
        cols = np.repeat(np.arange(4), 9)

        self.declare_partials(out_name, in_name, rows=rows, cols=cols)

    def compute(self, inputs, outputs):
        shape = self.options['shape']
        in_name = self.options['in_name']
        out_name = self.options['out_name']

        s = inputs[in_name][0][0]
        x = inputs[in_name][0][1]
        y = inputs[in_name][0][2]
        z = inputs[in_name][0][3]

        outputs[out_name] = np.array(
            [[1 - 2 * (y**2 + z**2), 2 * (x * y - s * z), 2 * (x * z + s * y)],
             [2 * (x * y + s * z), 1 - 2 * (x**2 + z**2), 2 * (y * z - s * x)],
             [2 * (x * z - s * y), 2 * (y * z + s * x),
              1 - 2 * (x**2 + y**2)]])

    def compute_partials(self, inputs, partials):
        in_name = self.options['in_name']
        out_name = self.options['out_name']

        s = inputs[in_name][0][0]
        x = inputs[in_name][0][1]
        y = inputs[in_name][0][2]
        z = inputs[in_name][0][3]

        partials[out_name, in_name] = np.array(
            [[0, -2 * z, 2 * y, 2 * z, 0, -2 * x, -2 * y, 2 * x, 0],
             [0, 2 * y, 2 * z, 2 * y, -4 * x, -2 * s, 2 * z, 2 * s, -4 * x],
             [-4 * y, 2 * x, 2 * s, 2 * x, 0, 2 * z, -2 * s, 2 * z, -4 * y],
             [-4 * z, -2 * s, 2 * x, 2 * s, -4 * z, 2 * y, 2 * x, 2 * y,
              0]]).flatten()


if __name__ == '__main__':
    from openmdao.api import Problem, IndepVarComp

    shape = (1, 4)

    prob = Problem()

    comp = IndepVarComp()
    comp.add_output('quat',
                    val=np.array([0, 0.267, 0.5345, 0.80178]).reshape(1, 4))
    prob.model.add_subsystem('ivc', comp, promotes=['*'])

    comp = Quat2RotComp(
        shape=shape,
        val=np.array([0, 0.267, 0.5345, 0.80178]).reshape(1, 4),
        in_name='quat',
        out_name='rot',
    )
    prob.model.add_subsystem('quat2rot_comp', comp, promotes=['*'])

    prob.setup()
    prob.run_model()
    prob.check_partials(compact_print=True)

    prob.model.list_inputs(print_arrays=True)
    prob.model.list_outputs(print_arrays=True)
