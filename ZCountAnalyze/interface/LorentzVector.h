#ifndef common_LorentzVector_h
#define common_LorentzVector_h

#include <TLorentzVector.h>

template<typename T>
TLorentzVector convert(T other) {
    return TLorentzVector(other.X(), other.Y(), other.Z(), other.T());
}

template<typename T>
TLorentzVector convert(T* other) {
    return TLorentzVector(other->X(), other->Y(), other->Z(), other->T());
}

template<typename T1, typename T2>
TLorentzVector sum(T1 first, T2 second) {
    return convert(first)+convert(second);
}

template<typename T1, typename T2>
double deltaR(T1 first, T2 second) {
    return convert(first).DeltaR(convert(second));
}

#endif // ifndef common_LorentzVector_h