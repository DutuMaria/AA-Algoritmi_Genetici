#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

int dimensiunePopulatie, limitaInferioara, limitaSuperioara, a, b, c, precizie, nrEtape, lungimeCromozom;
double probMutatie, probRecombinare;

ifstream fin ("input.in");
ofstream fout ("evolutie.out");

class Individ {
private:
    int valoareInt;
    double valoareDouble;
    double fitness;
    std::vector<int> cromozom;

public:
    Individ(int valoareInt, double valoareDouble, double fitness, const vector<int> &cromozom);

    int getValoareInt() const;
    void setValoareInt(int valInt);
    double getValoareDouble() const;
    void setValoareDouble(double valDouble);
    double getFitness() const;
    void setFitness(double fitness);
    const std::vector<int> &getCromozom() const;
    void setCromozom(const std::vector<int> &cromozomIndivid);

};


Individ::Individ(int valoareInt, double valoareDouble, double fitness, const vector<int> &cromozom) : valoareInt(
        valoareInt), valoareDouble(valoareDouble), fitness(fitness), cromozom(cromozom) {}

int Individ::getValoareInt() const { return valoareInt; }
void Individ::setValoareInt(int valInt) { Individ::valoareInt = valInt; }
double Individ::getValoareDouble() const { return valoareDouble; }
void Individ::setValoareDouble(double valDouble) { Individ::valoareDouble = valDouble; }
double Individ::getFitness() const { return fitness; }
void Individ::setFitness(double fitnessIndivid) { Individ::fitness = fitnessIndivid; }
const std::vector<int> &Individ::getCromozom() const { return cromozom; }
void Individ::setCromozom(const std::vector<int> &cromozomIndivid) { Individ::cromozom = cromozomIndivid; }

/*
 * -------------------------------------------------------------------------------
 */

// Functie care calculeaza 10^precizie
int calculareValTranslatie(){
    int valTranslatie = 1;
    int copiePrecizie = precizie;

    while(copiePrecizie > 0){
        valTranslatie *= 10;
        copiePrecizie--;
    }

    return valTranslatie;
}

int calculareLungimeCromozom() {
    int lungime = ceil(log2((b-a) * pow(10, precizie)));
    return lungime;
}

int generareValIntDinIntervalTranslatat() {
    int valTranslatie = calculareValTranslatie();

    /*
     * EXEMPLU
     * DACA avem - precizia = 10
     *           - interval [-1, 2]
     *           - limInf = -1, limSum = 2
     * => limInfTranslatata = -1 * 10 = -10
     * => limSumTranslatata = 2 * 10 = 20
     * scad limInfTranslatata din ambele (obs: prima data scad din limSup)
     * => limSupTranslatata = 20 - (-10) = 30
     * => limInfTranslatata = -10 - (-10) = 0
     */

    int limInfTranslatata = limitaInferioara * valTranslatie;
    int limSupTranslatata = limitaSuperioara * valTranslatie;

    limSupTranslatata -= limInfTranslatata;
    limInfTranslatata -= limInfTranslatata;

    std::random_device rd;
    std::mt19937 mt(rd());

    //trebuie minim pentru ca limSupTranslatata poate sa fie mai mare decat cel mai numar pe care il putem scrie pe lungimeCromozom biti
    limSupTranslatata = min(limSupTranslatata, (int)pow(2, lungimeCromozom)-1);
    std::uniform_int_distribution<> dist(limInfTranslatata, limSupTranslatata);

    int valIntRandom = dist(mt);
    return valIntRandom;
}

pair<int, double> generareValoariIndivid() {
    int valIntRandom = generareValIntDinIntervalTranslatat();
    int valTranslatie = calculareValTranslatie();
    double valoareIndivid = (double)(valIntRandom + limitaInferioara * valTranslatie) / valTranslatie;

    return {valIntRandom,valoareIndivid};
}

// fitness function => -X^2 + X + 2
double calculareFitnessIndivid(const double val) {
    double x = val;
    double fitnessIndivid = a * x * x + b * x + c;

    return fitnessIndivid;
}

double calculareFitnessTotal(const vector<Individ> &populatie){
    double fitnessTotal = 0;

    for(int i = 0; i < dimensiunePopulatie; i++){
        fitnessTotal += populatie[i].getFitness();
    }
    return fitnessTotal;
}

// calculare cromozom
std::vector<int> decimalToBinary(const int val) {
    int valIntRandom = val;
    vector<int> cromozomIndivid;
    int cnt = 0;
    while (valIntRandom > 0){
        cromozomIndivid.push_back(valIntRandom & 1);
        valIntRandom >>= 1;
        cnt ++;
    }

    int bitiRamasi = lungimeCromozom - cnt;

    while(bitiRamasi > 0){
        cromozomIndivid.push_back(0);
        bitiRamasi--;
    }
    reverse(cromozomIndivid.begin(), cromozomIndivid.end());

    return cromozomIndivid;
}

vector<Individ> generarePrimaPopulatie(){
    vector<Individ> populatie;
    for (int i = 0; i < dimensiunePopulatie; i++){
        pair<int, double> valoariIndivid = generareValoariIndivid();
        int valInt = valoariIndivid.first;
        double valDouble = valoariIndivid.second;
        double fitness = calculareFitnessIndivid(valDouble);
        vector<int> cromozom = decimalToBinary(valInt);
        Individ individ = Individ(valInt, valDouble, fitness, cromozom);
        populatie.push_back(individ);
    }

    return populatie;
}

Individ calculareIndividElitist(const vector<Individ> &populatie){
    double fitnessMax = populatie[0].getFitness();
    int cnt = 0;
    for (int i = 1; i < populatie.size(); i++){
        Individ individ = populatie[i];
        if (fitnessMax < individ.getFitness()){
            fitnessMax = individ.getFitness();
            cnt = i;
        }
    }

    Individ individElitist = populatie[cnt];

    return individElitist;
}

vector<double> calculareProbabilitatiSelectie(const vector<Individ> &populatie){
    double fitnessTotal = calculareFitnessTotal(populatie);
    vector<double> probabilitatiSelectie;
    for (const Individ& individ : populatie){
        double probSelectie = individ.getFitness() / fitnessTotal;
        probabilitatiSelectie.push_back(probSelectie);
    }

    return probabilitatiSelectie;
}

// suma cumulativa
vector<double> calculareIntervaleProbSelectie(const vector<double>& probabilitatiSelectie){
    vector<double> intervale;
    double sum = 0;
    intervale.push_back(0);
    for(double probabilitate : probabilitatiSelectie){
        sum += probabilitate;
        intervale.push_back(sum);
    }

    return intervale;
}

// cautare binara
int cautareInterval(const vector<double> &intervaleProbSel, const double &nr, int stanga, int dreapta){
    if (nr <= intervaleProbSel[stanga]){
        return stanga;
    }

    if (nr >= intervaleProbSel[dreapta]){
        return dreapta + 1;
    }

    if (stanga < dreapta){
        int mijl = (stanga + dreapta) / 2;

        if (intervaleProbSel[mijl] <= nr and intervaleProbSel[mijl + 1] > nr){
            return mijl + 1;
        } else if(intervaleProbSel[mijl] <= nr and intervaleProbSel[mijl +1] <= nr){
            return cautareInterval(intervaleProbSel, nr, mijl + 2, dreapta);
        } else {
            return cautareInterval(intervaleProbSel, nr, stanga, mijl - 1);
        }
    }
}

/*
 * => Procesul de selectie (putem alege un individ de mai multe ori)
 *      - selectam doar dimensiunePopulatie-1 indivizi pentru ca adaug la final elitistul
 *      - generez un nunar intre [0, nrIntervale]  unde nrIntervale = intervalePobSel.size() - 1
 *      - caut pozitia numarului in interval
 *      - selectez individul de pe pozitia  poz-1 (exemplu pt obs: ca individul 1 e pe pozitia 0!!)
 */
vector<Individ> selecteazaIndivizi(const int &etapa, const vector<Individ> &populatie, const vector<double> &intervalePobSel){
    int indiviziDeSelectat = dimensiunePopulatie - 1;
    vector<Individ> indiviziSelectati;
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0, 1);

    while (indiviziDeSelectat > 0){
        double nr = dist(mt);
        int poz = cautareInterval(intervalePobSel, nr, 0, (int)intervalePobSel.size() - 1);

        if (etapa == 0){
            fout << "\tu= " << nr << " selectam cromozomul " << poz << "\n";
        }
        indiviziSelectati.push_back(populatie[poz - 1]); // exemplu: individul 1 e pe pozitia 0!!
        indiviziDeSelectat--;
    }

    return indiviziSelectati;
}

void afisareCromozom(const vector<int> &cromozom){
    for (auto gena : cromozom){
        fout << gena;
    }
}

/*
 * Functia returneaza un vector cu pozitiile indivizilor care trebuie recombinati
 *      - generez un numar aleator pentru fiecare individ
 *      - daca numarul generat < probabilitatea de recombinare
 *          - adaug pozitia individului la vectorul cu pozitiile indivizilor care trebuie recombinati
 */
vector<int> selecteazaIndiviziDeIncrucisat(const int &etapa, const vector<Individ> &indiviziSelectati){
    vector<int> indiviziSelectatiIncrucisare;
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0, 1);

    if (etapa == 0){ fout << "\nProbabilitatea de incrucisare: " << probRecombinare << "\n"; }
    for (int i = 0; i < indiviziSelectati.size(); i++){
        Individ individ = indiviziSelectati[i];
        vector<int> cromozom = individ.getCromozom();
        double nr = dist(mt);
        if (nr < probRecombinare){
            indiviziSelectatiIncrucisare.push_back(i);
            if (etapa == 0){
                fout << "\t" << i + 1 << ": ";
                afisareCromozom(cromozom);
                fout << " u= " << nr << " < " << probRecombinare << " participa\n";
            }
        } else{
            if (etapa == 0){
                fout << "\t" << i + 1 << ": ";
                afisareCromozom(cromozom);
                fout << " u= " << nr << " \n";
            }
        }
    }
    if (etapa == 0){ fout << " \n"; }

    return indiviziSelectatiIncrucisare;
}

int binaryToDecimal(vector<int> cromozom){
    reverse(cromozom.begin(), cromozom.end());
    int valInt = 0;

    for (int i = 0; i < lungimeCromozom; i++){
        valInt += (cromozom[i] & 1) * (1 << i);
    }

    return valInt;
}

double intToDouble(const int &valInt){
    int valTranslatie = calculareValTranslatie();
    double valDouble = (double)(valInt + limitaInferioara * valTranslatie) / valTranslatie;

    return valDouble;
}


/*
 * Procesul de recombinare
 *  => returnez un vector de perechi (pozitieIndivid, Individ) cu indivizii rezultati in urma recombinarii
 */
vector<pair<int,Individ>> recombinare(const int &etapa, const vector<Individ> & indiviziSelectati, vector<int> &indiviziDeIncrucisat){
    int nrIndiviziDeIncrucisat = indiviziDeIncrucisat.size();
    vector<pair<int, Individ>> indiviziRecombinati;

    shuffle(indiviziDeIncrucisat.begin(), indiviziDeIncrucisat.end(), std::mt19937(std::random_device()()));

    if (nrIndiviziDeIncrucisat % 2 == 1){
        nrIndiviziDeIncrucisat --;
        int indiceIndividNeschimbat = indiviziDeIncrucisat[nrIndiviziDeIncrucisat];
        indiviziRecombinati.push_back({indiceIndividNeschimbat, indiviziSelectati[indiceIndividNeschimbat]});
    }

    for (int i = 0; i < nrIndiviziDeIncrucisat - 1; i += 2){
        int indiceIndivid1 = indiviziDeIncrucisat[i], indiceIndivid2 = indiviziDeIncrucisat[i + 1];
        Individ individ1 = indiviziSelectati[indiceIndivid1], individ2 = indiviziSelectati[indiceIndivid2];
        vector<int> cromozom1 = individ1.getCromozom(), cromozom2 = individ2.getCromozom();
        int punctRupere;
        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_int_distribution<> dist(0, lungimeCromozom);
        punctRupere = dist(mt);

        if (etapa == 0){
            fout << "Recombinare dintre cromozomul " << indiceIndivid1 + 1 << " cu cromozomul " << indiceIndivid2 + 1<< ":\n";
            afisareCromozom(cromozom1);
            fout << " ";
            afisareCromozom(cromozom2);
            fout << " punct " << punctRupere << "\n";
            fout << "Rezultat\t";
        }

        vector<int> cromozomRecombinat1(lungimeCromozom), cromozomRecombinat2(lungimeCromozom);

        for (int j = 0; j < lungimeCromozom; j++){
            if (j < punctRupere){
                cromozomRecombinat1[j] = cromozom2[j];
                cromozomRecombinat2[j] = cromozom1[j];
            } else {
                cromozomRecombinat1[j] = cromozom1[j];
                cromozomRecombinat2[j] = cromozom2[j];
            }
        }

        if (etapa == 0){
            afisareCromozom(cromozomRecombinat1);
            fout << " ";
            afisareCromozom(cromozomRecombinat2);
            fout << "\n\n";
        }

        int valInt1 = binaryToDecimal(cromozomRecombinat1);
        double valDouble1 = intToDouble(valInt1);
        double fitness1 = calculareFitnessIndivid(valDouble1);
        Individ individRecombinat1 = Individ(valInt1, valDouble1, fitness1, cromozomRecombinat1);
        indiviziRecombinati.push_back({indiceIndivid1, individRecombinat1});

        int valInt2 = binaryToDecimal(cromozomRecombinat2);
        double valDouble2 = intToDouble(valInt2);
        double fitness2 = calculareFitnessIndivid(valDouble2);
        Individ individRecombinat2 = Individ(valInt2, valDouble2, fitness2, cromozomRecombinat2);
        indiviziRecombinati.push_back({indiceIndivid2, individRecombinat2});
    }

    return indiviziRecombinati;
}

/*
 * Procesul de mutatie
 *  => returnez un vector de perechi (pozitieIndivid, Individ) cu indivizii rezultati in urma mutatiei
 */
vector<pair<int, Individ>> mutatie(const vector<Individ> &populatieRecombinata){
    vector<pair<int, Individ>> indiviziModificati;
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0, 1);

    for (int i = 0; i < dimensiunePopulatie - 1; i++){
        bool esteModificat = false;
        Individ individ = populatieRecombinata[i];
        vector<int> cromozom = individ.getCromozom();
        vector<int> cromozomNou;

        for (int j = 0; j < lungimeCromozom; j++){
            double nr = dist(mt);

            if(nr < probMutatie){
                esteModificat = true;
                cromozomNou.push_back(cromozom[j] ^ 1); // 0 ^ 1 = 1 si 1 ^ 1 = 0
            } else {
                cromozomNou.push_back(cromozom[j]);
            }

        }

        if(esteModificat){
            int valInt = binaryToDecimal(cromozomNou);
            double valDouble = intToDouble(valInt);
            double fitness = calculareFitnessIndivid(valDouble);
            Individ individModificat = Individ(valInt, valDouble, fitness, cromozomNou);
            indiviziModificati.push_back({i, individModificat});
        }
    }

    return indiviziModificati;
}

void afisareIndivid(const int &indice, const Individ &individ){
    fout << "\t" << indice + 1 << ": ";
    const vector<int>& cromozom = individ.getCromozom();
    afisareCromozom(cromozom);
    fout << " x= " << individ.getValoareDouble() << " f= " << individ.getFitness() << "\n";
}

void afisarePopulatieInitiala(const vector<Individ> &primaPopulatie){
    fout << "Populatia initiala\n";
    for (int j = 0; j < dimensiunePopulatie; j++) {
        Individ individ = primaPopulatie[j];
        afisareIndivid(j, individ);
    }
}

void afisareProbabilitatiSelectie(const vector<double> &probabilitatiSelectie, const vector<Individ> &populatie){
    fout << "\nProbabilitati selectie\n";
    for (int j = 0; j < dimensiunePopulatie; j++) {
        Individ individ = populatie[j];
        fout << "\tcromozom\t" << j + 1 << " probabilitate " << probabilitatiSelectie[j] << "\n";
    }
}

void afisareIntervaleProbSel(const vector<double> &intervaleProbSelectie){
    fout << "\nIntervale probabilitati selectie\n";
    for (auto nr : intervaleProbSelectie){
        fout << nr << " ";
    }
    fout << "\n\n";
}

void afisareDupaSelectie(vector<Individ> & indiviziSelectati){
    fout << "\nDupa selectie\n";
    for (int j = 0; j < dimensiunePopulatie - 1; j++){
        Individ individ = indiviziSelectati[j];
        afisareIndivid(j, individ);
    }
}

void afisareDupaRecombinare(const vector<Individ> &populatieDupaRecombinare){
    fout << "Dupa recombinare:\n";
    for (int j = 0; j < dimensiunePopulatie - 1; j++) {
        Individ individ = populatieDupaRecombinare[j];
        afisareIndivid(j, individ);
    }
}

void afisareCromozomiModifDupaMutatie(const vector<pair<int, Individ>> &indiviziModificati){
    fout << "\nProbabilitate de mutatie pentru fiecare gena " << probMutatie << "\n";
    fout << "Au fost modificati cromozomii:\n";

    for(const auto& pereche : indiviziModificati){
        fout << "\t" << pereche.first + 1 << "\n";
    }
}

void afisareDupaMutatie(const vector<Individ> &populatieDupaMutatie){
    fout << "\nDupa mutatie:\n";
    for (int j = 0; j < dimensiunePopulatie - 1; j++) {
        Individ individ = populatieDupaMutatie[j];
        afisareIndivid(j, individ);
    }
}

void afisareEvolutieMaxim(const vector<double> &evolutieFitnessMaxim){
    fout << "\nEvolutia maximului\n";
    for (double fitness : evolutieFitnessMaxim){
        fout << setprecision(16) << fitness << "\n";
    }
}

void afisareValoareMediePerformanta(const vector<double> &valoriMediePerformanta){
    fout << "\nValoare medie performanta\n";
    for (double performanta : valoriMediePerformanta){
        fout << setprecision(16) <<  performanta << "\n";
    }
}


int main() {
    fin >> dimensiunePopulatie >> limitaInferioara >> limitaSuperioara >> a >> b >> c >> precizie >> probRecombinare >> probMutatie >> nrEtape;
    lungimeCromozom = calculareLungimeCromozom();
    fout << "Evolutie.out\n\n";
    vector<Individ> populatie = generarePrimaPopulatie();
    vector<double> evolutieMaxim;
    vector<double> valoriMediePerformanta;

    for (int i = 0; i < nrEtape; i++){
        if (i == 0) { afisarePopulatieInitiala(populatie); }

        vector<double> probabilitatiSelectie = calculareProbabilitatiSelectie(populatie);
        if (i == 0) { afisareProbabilitatiSelectie(probabilitatiSelectie, populatie); }

        vector<double> intervaleProbSelectie = calculareIntervaleProbSelectie(probabilitatiSelectie);
        if (i == 0) { afisareIntervaleProbSel(intervaleProbSelectie); }

        vector<Individ> indiviziSelectati = selecteazaIndivizi(i, populatie, intervaleProbSelectie);
        if (i == 0) { afisareDupaSelectie(indiviziSelectati); }

        vector<int> indiviziDeIncrucisat = selecteazaIndiviziDeIncrucisat(i, indiviziSelectati);
        vector<pair<int, Individ>> indiviziRecombinati = recombinare(i, indiviziSelectati, indiviziDeIncrucisat);
        vector<Individ> populatieDupaRecombinare = indiviziSelectati;

        for (int j = 0; j < indiviziRecombinati.size(); j++){
            int indiceIndividRecombinat = indiviziRecombinati[j].first;
            Individ individRecombinat = indiviziRecombinati[j].second;
            populatieDupaRecombinare[indiceIndividRecombinat] = individRecombinat;
        }

        if (i == 0) { afisareDupaRecombinare(populatieDupaRecombinare); }

        vector<pair<int, Individ>> indiviziModificatiDupaMutatie = mutatie(populatieDupaRecombinare);
        if (i == 0) { afisareCromozomiModifDupaMutatie(indiviziModificatiDupaMutatie); }
        vector<Individ> populatieDupaMutatie = populatieDupaRecombinare;

        for (int j = 0; j < indiviziModificatiDupaMutatie.size(); j++){
            int indiceIndividModificat = indiviziModificatiDupaMutatie[j].first;
            Individ individModificat = indiviziModificatiDupaMutatie[j].second;
            populatieDupaMutatie[indiceIndividModificat] = individModificat;
        }

        if (i == 0) { afisareDupaMutatie( populatieDupaMutatie); }

        vector<Individ> populatieNoua = populatieDupaMutatie;
        Individ individElitist = calculareIndividElitist(populatie);
        populatieNoua.push_back(individElitist);
        populatie = populatieNoua;
        evolutieMaxim.push_back(individElitist.getFitness());
        valoriMediePerformanta.push_back(calculareFitnessTotal(populatie) / dimensiunePopulatie);
    }

    afisareEvolutieMaxim(evolutieMaxim);
    afisareValoareMediePerformanta(valoriMediePerformanta);

    return 0;
}
