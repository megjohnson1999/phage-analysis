from abc import abstractmethod
from enum import Enum
from typing import FrozenSet, List, Set, TypeVar

TSettingsEnumBase = TypeVar("TSettingsEnumBase", bound="SettingsEnumBase")


class SettingsEnumBase(Enum):
    @classmethod
    def choices(cls) -> List[str]:
        return sorted(item.item_to_choice() for item in cls)

    @classmethod
    def all(cls) -> FrozenSet[TSettingsEnumBase]:
        skip = cls.skip_for_all()
        return frozenset(item for item in cls if item not in skip)

    @classmethod
    def parse_choices_list(cls, choices: List[str]) -> List[TSettingsEnumBase]:
        return cls._parse_choices_into(choices, list)

    @classmethod
    def parse_choices_set(cls, choices: List[str]) -> Set[TSettingsEnumBase]:
        return cls._parse_choices_into(choices, set)

    @classmethod
    def _parse_choices_into(cls, choices: str, container) -> List[TSettingsEnumBase]:
        return container(cls.parse_choice(choice) for choice in choices)

    @classmethod
    def parse_choice(cls, choice: str) -> TSettingsEnumBase:
        return cls[choice.replace("-", "_").upper()]

    @classmethod
    @abstractmethod
    def skip_for_all(cls) -> FrozenSet[TSettingsEnumBase]:
        return frozenset()

    def item_to_choice(self) -> str:
        return self.name.replace("_", "-").lower()

    def __str__(self):
        return self.item_to_choice()
