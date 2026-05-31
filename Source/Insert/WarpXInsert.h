#pragma once

namespace Insert
{

void Initialize ();
void BeforeStep ();
void ParticleInjection ();
void SetBoundaryPhi ();
void SetPhiGuards ();
void BeforeCollision (int step);
void AfterCollision (int step);
void AfterDiagnostics ();

} // namespace Insert
