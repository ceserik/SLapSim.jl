function tire(longitudnalForce, slipAngle,optmizationProblem,Fz)
    lateralForce = slipAngle*Fz*constant

    optmizationProblem.subjectTo(lateralForce^2 + longitudnalForce^2 <= 50000)
    optmizationProblem.subjectTo(slipAngle<5)
    
    return [lateralForce,optimizationProblem]
end