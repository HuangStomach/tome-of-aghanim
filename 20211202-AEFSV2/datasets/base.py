from abc import ABCMeta
from abc import abstractmethod

class Base:
    __metaclass__= ABCMeta

    @abstractmethod
    def drugs(self):
        pass
