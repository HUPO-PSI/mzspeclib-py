"""Tools for interacting with controlled vocabularies and biomedical ontologies."""

import logging
from typing import Dict

from psims.controlled_vocabulary import Entity, ControlledVocabulary
from psims.controlled_vocabulary.controlled_vocabulary import load_uo, load_unimod, load_psims

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class _VocabularyResolverMixin(object):
    """A base class for looking up terms in a set of controlled vocabularies."""

    default_cv_loader_map = {
        "MS": load_psims,
        "UO": load_uo,
        "UNIMOD": load_unimod,
    }

    controlled_vocabularies: Dict[str, ControlledVocabulary]

    def __init__(self, *args, **kwargs):
        self.controlled_vocabularies = dict()
        super().__init__(*args, **kwargs)

    def load_cv(self, name: str) -> ControlledVocabulary:
        """
        Retrieve a controlled vocabulary by name.

        This may query the internet or on-disk cache when called for the
        first time for each name, and the result is cached in memory.
        There may be inconsistencies from one source to the next.

        Parameters
        ----------
        name : str
            The name to load

        Returns
        -------
        :class:`psims.controlled_vocabular.ControlledVocabulary`
        """
        if name in self.controlled_vocabularies:
            return self.controlled_vocabularies[name]
        self.controlled_vocabularies[name] = self.default_cv_loader_map[name]()
        return self.controlled_vocabularies[name]

    def find_term_for(self, curie: str) -> Entity:
        try:
            name, _id = curie.split(":", 1)
        except ValueError as err:
            raise ValueError(f"Could not parse {curie} into source and accession. Is it a CURIE?") from None
        cv = self.load_cv(name)
        term = cv[curie]
        return term

    def find_term_by_name(self, name: str) -> Entity:
        cv_keys = ['MS', 'UNIMOD']
        for cv_key in cv_keys:
            cv = self.load_cv(cv_key)
            try:
                term = cv[name]
                return term
            except KeyError:
                continue
        raise KeyError(name)


class ControlledVocabularyResolver(_VocabularyResolverMixin):

    def name_to_curie(self, name: str) -> str:
        term = self.find_term_by_name(name)
        return term.id

    def attribute_syntax(self, name: str) -> str:
        """Given an identifier, format it for textual representation"""
        if self.is_curie(name):
            if "|" in name:
                return name
            term = self.find_term_for(name)
        else:
            term = self.find_term_by_name(name)
        return f"{term.id}|{term.name}"

    def is_curie(self, name: str) -> bool:
        """Test if a string is a CURIE"""
        if ':' not in name:
            return False
        prefix, _rest = name.split(":")
        if prefix in self.default_cv_loader_map:
            return True
        return False
