from abc import ABC, abstractmethod
from pathlib import Path


class RuleInterface(ABC):
    @property
    @abstractmethod
    def name(self) -> str: ...

    @property
    @abstractmethod
    def lineno(self) -> int: ...

    @property
    @abstractmethod
    def snakefile(self) -> Path: ...
