"""Tested *material-point* ductile-damage prototype; NOT coupled to CalculiX."""
from dataclasses import dataclass

@dataclass
class DuctilePoint:
    accumulated_plastic_strain: float = 0.0
    initiation_index: float = 0.0
    damage: float = 0.0
    plastic_displacement_after_initiation: float = 0.0
    failed: bool = False

    def update(self, delta_peeq, triaxiality, characteristic_length, fracture_curve, failure_displacement):
        """Accumulate strain-weighted initiation and linear post-initiation displacement.

        fracture_curve: sorted (triaxiality, fracture_plastic_strain) pairs.
        This is a simplified material-point algorithm, not Abaqus' integrated damage law.
        """
        import numpy as np
        if delta_peeq < 0 or characteristic_length <= 0 or failure_displacement <= 0:
            raise ValueError('Invalid damage increment/length')
        if self.failed: return self
        etas, eps = zip(*sorted(fracture_curve))
        critical=float(np.interp(triaxiality,etas,eps))
        if critical <= 0: raise ValueError('Fracture strain must be positive')
        old=self.initiation_index
        inc=delta_peeq/critical
        self.initiation_index += inc
        self.accumulated_plastic_strain += delta_peeq
        # Only the increment after crossing initiation contributes to softening.
        post_fraction = 0 if inc==0 else min(1.0,max(0.0,(old+inc-1.0)/inc))
        self.plastic_displacement_after_initiation += characteristic_length*delta_peeq*post_fraction
        self.damage=min(1.0,self.plastic_displacement_after_initiation/failure_displacement)
        self.failed=self.damage>=1.0
        return self
