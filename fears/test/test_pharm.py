import pytest
from fears.utils.pharm import *
from fears.population import Population, PopParams

@pytest.fixture
def pop():
    pop = Population(fitness_data='random')
    return pop

def test_pharm_eqn(pop):
    p = pharm_eqn(pop, t = 100)
    assert isinstance(p, float)

def test_pharm_eqn2(pop):
    p = pharm_eqn(pop, t = 100)
    assert p >=0

def test_convolve_pharm(pop):
    u = [0,1,0,1,0,0,0,0,0,0,0,0,1,1,1,1,0,0]
    p = convolve_pharm(pop,u)
    assert len(p) == pop.n_timestep

def test_convolve_pharm2(pop):
    u = [0,1,0,1,0,0,0,0,0,0,0,0,1,1,1,1,0,0]
    p = convolve_pharm(pop,u)
    assert isinstance(p[0], float)

def test_convolve_pharm3(pop):
    u = [0,1,0,1,0,0,0,0,0,0,0,0,1,1,1,1,0,0]
    p = convolve_pharm(pop,u)
    bools = [i >=0 for i in p]
    assert all(bools)

def test_get_impulses(pop):
    i = gen_impulses(pop)
    assert len(i) == pop.n_timestep

def test_get_impulses2(pop):
    imp = gen_impulses(pop)
    bools = [i ==0 or i ==1 for i in imp]
    assert all(bools)

#def test_gen_on_off(pop):
#    v = gen_on_off_regimen(pop)

def test_gen_curves(pop):
    c = gen_curves(pop)
    curve = c[0]
    assert len(curve) == pop.n_timestep

def test_gen_passage_drug(pop):
    p = gen_passage_drug_protocol(pop)
    assert len(p) == pop.n_timestep