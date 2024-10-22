import warnings

from typing import Iterable, KeysView, ItemsView, Optional, ValuesView, Dict, MutableMapping, Mapping

from pyteomics import proforma

from mzspeclib.attributes import AttributedEntity, IdentifiedAttributeManager, AttributeManagedProperty, AttributeProxy, AttributeGroupFacet
from mzspeclib.const import (ANALYTE_MIXTURE, CHARGE_STATE, PROFORMA_ION, PROFORMA_SEQ, STRIPPED_PEPTIDE_SEQ, FIRST_ANALYTE_KEY, FIRST_INTERPRETATION_KEY)


class _AnalyteMappingProxy(Mapping[str, 'Analyte']):
    parent: Mapping[str, 'Interpretation']

    def __init__(self, parent):
        self.parent = parent

    def __getitem__(self, key):
        for group_id, group in self.parent.items():
            if key in group:
                return group[key]
        raise KeyError(key)

    def __iter__(self):
        keys = set()
        for group_id, group in self.parent.items():
            k_ = set(group.keys())
            keys.update(k_)
        return iter(keys)

    def __len__(self):
        return sum([len(v) for v in self.parent.values()])

    def __repr__(self):
        d = dict(self)
        return f"{self.__class__.__name__}({d})"


class InterpretationCollection(MutableMapping[str, 'Interpretation']):
    """
    A mutable mapping for :class:`Interpretation`s that also exposes a shared pool of
    :class:`Analyte` members.
    """

    __slots__ = ('interpretations', )

    interpretations: Dict[str, 'Interpretation']

    def __init__(self, interpretations=None):
        if interpretations is None:
            interpretations = {}
        self.interpretations = interpretations

    def get_interpretation(self, interpretation_id) -> 'Interpretation':
        """Get the interpretation given by the key ``interpretation_id``"""
        return self.interpretations[str(interpretation_id)]

    def add_interpretation(self, interpretation: 'Interpretation'):
        """
        Add an :class:`Interpretation` to the collection.

        The key used will be :attr:`Interpretation.id`.
        """
        self.set_interpretation(str(interpretation.id), interpretation)

    def set_interpretation(self, key, interpretation: 'Interpretation'):
        """Set the interpretation for ``key`` to ``interpretation``"""
        self.interpretations[str(key)] = interpretation

    def __getitem__(self, key) -> 'Interpretation':
        return self.get_interpretation(key)

    def __setitem__(self, key, value: 'Interpretation'):
        self.set_interpretation(key, value)

    def __delitem__(self, key):
        del self.interpretations[str(key)]

    def __len__(self):
        return len(self.interpretations)

    def __contains__(self, key):
        return key in self.interpretations

    def __iter__(self):
        return iter(self.interpretations)

    def keys(self) -> KeysView[str]:
        return self.interpretations.keys()

    def values(self) -> ValuesView['Interpretation']:
        return self.interpretations.values()

    def items(self) -> ItemsView[str, 'Interpretation']:
        return self.interpretations.items()

    @property
    def analytes(self) -> '_AnalyteMappingProxy':
        """A facade that exposes a :class:`Mapping[str, Analyte]` interface"""
        return _AnalyteMappingProxy(self)

    def __repr__(self):
        d = dict(self)
        return f"{self.__class__.__name__}({d})"


class Interpretation(AttributedEntity, MutableMapping):
    """
    An interpretation of a :class:`~.Spectrum` with one or more :class:`~.Analyte`
    members.

    Attributes
    ----------
    id : str
        The identifier for the interpretation
    analytes : dict[str, :class:`Analyte`]
        The analytes which are part of interpretations.
    member_interpretations : dict[str, :class:`InterpretationMember`]
        The interpretation details which are associated with specific :class:`Analyte` members.
    """

    __slots__ = ('id', 'analytes', 'member_interpretations')

    id: str
    analytes: Dict[str, 'Analyte']
    member_interpretations: Dict[str, 'InterpretationMember']

    def __init__(self, id, attributes: Iterable = None, analytes: Dict = None, member_interpretations: Dict = None):
        self.id = str(id)
        self.analytes = analytes or {}
        self.member_interpretations = member_interpretations or {}
        super(Interpretation, self).__init__(attributes)

    def _update_mixture_members_term(self):
        value = sorted(map(int, self.analytes.keys()))
        self.replace_attribute(ANALYTE_MIXTURE, value)

    def get_analyte(self, analyte_id) -> 'Analyte':
        """Retrieve an analyte by its identifier"""
        return self.analytes[str(analyte_id)]

    def add_analyte(self, analyte: 'Analyte'):
        """Add an analyte to the interpretation"""
        self.set_analyte(analyte.id, analyte)
        self._update_mixture_members_term()

    def set_analyte(self, key, analyte: 'Analyte'):
        """Set the analyte for ``key`` to ``analyte``"""
        key = str(key)
        self.analytes[key] = analyte
        self._update_mixture_members_term()

    def remove_analyte(self, analyte_id):
        """Remove the analyte for ``analyte_id``"""
        del self.analytes[str(analyte_id)]
        self._update_mixture_members_term()

    def has_analyte(self, analyte_id) -> bool:
        """Check if this interpretation includes ``analyte_id``"""
        return str(analyte_id) in self.analytes

    def get_member_interpretation(self, member_id) -> 'InterpretationMember':
        """Retrieve the :class:`InterpretationMember` for the ``member_id``"""
        return self.member_interpretations[str(member_id)]

    def add_member_interpretation(self, interpretation_member: 'InterpretationMember'):
        """Add an :class:`InterpretationMember` to the interpretation"""
        self._set_member_interpretation(interpretation_member.id, interpretation_member)

    def _set_member_interpretation(self, key, interpretation_member: 'InterpretationMember'):
        self.member_interpretations[str(key)] = interpretation_member

    def remove_member_interpretation(self, member_id):
        """Remove the :class:`InterpretationMember` for ``member_id``"""
        del self.member_interpretations[str(member_id)]

    def validate(self) -> bool:
        """
        Perform validation on each component to confirm this object is well formed.

        Returns
        -------
        bool
        """
        analyte_ids = set(self.analytes)
        member_ids = set(self.member_interpretations)
        valid = True
        if not (analyte_ids >= member_ids):
            warnings.warn(
                f"Interpretation has InterpretationMembers {member_ids - analyte_ids} lacking Analytes")
            valid = False
        return valid

    def __getitem__(self, key) -> 'Analyte':
        return self.get_analyte(key)

    def __setitem__(self, key, value: 'Analyte'):
        self.set_analyte(key, value)

    def __delitem__(self, key):
        self.remove_analyte(key)

    def __iter__(self):
        return iter(self.analytes)

    def __len__(self):
        return len(self.analytes)

    def __repr__(self):
        d = dict(self)
        a = ''
        if self.attributes:
            a = f', {self.attributes}'
        return f"{self.__class__.__name__}({d}{a})"


class InterpretationMember(IdentifiedAttributeManager):
    """
    A collection of attributes associated with a specific :class:`Analyte`
    contained in an :class:`Interpretation`
    """

    __slots__ = ()


class ProteinDescription(AttributeProxy):
    """
    A protein associated with a peptide :class:`Analyte`.

    This is a read-only proxy.
    """

    accession = AttributeManagedProperty("MS:1000885|protein accession")
    name = AttributeManagedProperty("MS:1000886|protein name")
    missed_cleavages = AttributeManagedProperty[int]("MS:1003044|number of missed cleavages")
    cleavage_agent = AttributeManagedProperty("MS:1001045|cleavage agent name")
    number_of_enzymatic_termini = AttributeManagedProperty[int]("MS:1003048|number of enzymatic termini")
    flanking_n_terminal_residue = AttributeManagedProperty("MS:1001112|n-terminal flanking residue")
    flanking_c_terminal_residue = AttributeManagedProperty("MS:1001113|c-terminal flanking residue")


class Analyte(IdentifiedAttributeManager):
    """A molecule that is associated with an :class:`Interpretation`"""

    __slots__ = ()

    mass = AttributeManagedProperty[float]("MS:1001117|theoretical mass")
    proteins = AttributeGroupFacet[ProteinDescription](ProteinDescription)

    @property
    def peptide(self) -> Optional[proforma.ProForma]:
        """
        Read out the peptide sequence of the analyte if it is present.

        This probes the following attributes, in order:
            1. MS:1003270|proforma peptidoform ion notation
            2. MS:1000889|proforma peptidoform sequence
            3. MS:1000888|stripped peptide sequence

        Returns
        -------
        :class:`~pyteomics.proforma.ProForma` or :const:`None`
        """
        if self.has_attribute(PROFORMA_SEQ):
            peptide_seq = self.get_attribute(PROFORMA_SEQ)
            return proforma.ProForma.parse(peptide_seq)
        if self.has_attribute(PROFORMA_ION):
            peptide_seq = self.get_attribute(PROFORMA_ION)
            return proforma.ProForma.parse(peptide_seq)
        if self.has_attribute(STRIPPED_PEPTIDE_SEQ):
            peptide_seq = self.get_attribute(STRIPPED_PEPTIDE_SEQ)
            return proforma.ProForma.parse(peptide_seq)
        return None

    @property
    def charge(self) -> Optional[int]:
        """
        Read the analyte's charge state, if it is present.

        This probes the following attributes in order:
            1. MS:1000041|charge state
            2. MS:1003270|proforma peptidoform ion notation

        Returns
        -------
        :class:`int` or :const:`None`
        """
        if self.has_attribute(CHARGE_STATE):
            return self.get_attribute(CHARGE_STATE)
        elif self.has_attribute(PROFORMA_ION):
            ion_val = self.get_attribute(PROFORMA_ION)
            val = proforma.ProForma.parse(ion_val)
            return val.charge_state
        else:
            return None

    @charge.setter
    def charge(self, value):
        if value is not None:
            if self.has_attribute(CHARGE_STATE):
                self.replace_attribute(CHARGE_STATE, value)
            else:
                self.add_attribute(CHARGE_STATE, value)
        else:
            self.remove_attribute(CHARGE_STATE)