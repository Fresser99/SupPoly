from enum import Enum


class CompType(Enum):

    conventional=0,
    polymer=1


class Component:

    def __init__(self, formular, cas,type):
        self.Formular = formular
        self.CAS = cas
        self.type=type
