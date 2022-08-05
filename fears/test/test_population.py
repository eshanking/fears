from fears.population import Population, PopParams
import pytest

@pytest.fixture
def default_pop():
    p = Population()
    return p

@pytest.fixture
def pop_estimated_data():
    p = Population(fitness_data='estimate')
    return p

@pytest.fixture
def pop_random_data():
    p = Population(fitness_data='random',n_allele=2)
    return p

def test_default_data(default_pop):
    assert default_pop.n_allele == 4
    assert default_pop.n_genotype == 16
    assert len(default_pop.ic50) == 16
    assert len(default_pop.drugless_rates) == 16

def test_estimated_data(pop_estimated_data):
    assert pop_estimated_data.n_allele == 4
    assert pop_estimated_data.n_genotype == 16
    assert len(pop_estimated_data.ic50) == 16
    assert len(pop_estimated_data.drugless_rates) == 16

def test_estimated_data(pop_random_data):
    assert pop_random_data.n_allele == 2
    assert pop_random_data.n_genotype == 4
    assert len(pop_random_data.ic50) == 4
    assert len(pop_random_data.drugless_rates) == 4