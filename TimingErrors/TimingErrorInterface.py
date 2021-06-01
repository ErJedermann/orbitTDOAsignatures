import abc


class TimingError(abc.ABC):
    @abc.abstractmethod
    def __init__(self, mean, std):
        pass

    @abc.abstractmethod
    def __str__(self):
        pass

    @abc.abstractmethod
    def get_error(self, amount):
        pass

    @abc.abstractmethod
    def get_mean_std(self):
        pass
