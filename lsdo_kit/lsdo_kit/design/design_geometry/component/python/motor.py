from lsdo_kit.components.component import Component

class Motor(Component):
    def __init__(self, stp_name: str , alias: str ):

        super(Motor, self).__init__(stp_name, alias)
