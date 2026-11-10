#include <SPI.h>
#include <arduinoFFT.h>

// =====================================================
// CONFIGURAÇÃO
// =====================================================

#define SAMPLES 256
#define SAMPLING_FREQUENCY 3840.0

// MCP3208
#define MCP3208_CS 53

#define MCP_VREF 5.0
#define MCP_MAX 4095.0

// SPI
#define SPI_CLOCK 1000000UL   // 1 MHz

// =====================================================
// CANAIS DO MCP3208
// =====================================================

#define V1_CH 0
#define V2_CH 1
#define V3_CH 2

#define I1_CH 3
#define I2_CH 4
#define I3_CH 5

// =====================================================
// BUFFERS ADC
// =====================================================

uint16_t adc[6][SAMPLES];

// =====================================================
// FFT
// =====================================================

float vReal[SAMPLES];
float vImag[SAMPLES];

ArduinoFFT<float> FFT;

// =====================================================
// HARMÔNICAS DE PICO
// =====================================================

int harmPicoV[3];
float harmPicoAmpV[3];
float harmPicoFreqV[3];

int harmPicoI[3];
float harmPicoAmpI[3];
float harmPicoFreqI[3];

// =====================================================
// ESPECTRO DE TENSÃO
// =====================================================

float harmAmpV[3][10];
float harmFreqV[3][10];
int harmBinV[3][10];

// =====================================================
// ESPECTRO DE CORRENTE
// =====================================================

float harmAmpI[3][10];
float harmFreqI[3][10];
int harmBinI[3][10];

// =====================================================
// GANHOS
// =====================================================

float GainV = 505.1167;
float GainI = 571.404;

// =====================================================
// RESULTADOS
// =====================================================

float Vrms[3];
float Irms[3];

float PotAtiva[3];
float PotAparente[3];
float FP[3];

float Ptotal = 0;
float Stotal = 0;
float FPtotal = 0;

// =====================================================
// FFT AUXILIARES
// =====================================================

float freqFundamental = 0;
float faseFund[6];
float phi[3];
float DPF[3];

int binFundamental = 1;

// =====================================================
// HARMÔNICAS E THD
// =====================================================

float THDv[3];
float THDi[3];

float VrmsFund[3];
float IrmsFund[3];

float VrmsTotalFFT[3];
float IrmsTotalFFT[3];

float DF[3];

// =====================================================
// TEMPO
// =====================================================

unsigned long Ts_us;

// =====================================================
// LEITURA DO MCP3208
// =====================================================
//
// MCP3208 - leitura single-ended
//
// Byte 1:
// 00000110 | canal >> 2
//
// Byte 2:
// canal restante nos bits superiores
//
// Byte 3:
// recebe os bits restantes
//
// Retorno: 0 a 4095
// =====================================================

uint16_t lerMCP3208(uint8_t canal)
{
    uint8_t comando;
    uint8_t byte2;
    uint8_t byte3;

    canal &= 0x07;

    comando = 0x06 | (canal >> 2);
    byte2 = (canal & 0x03) << 6;

    digitalWrite(MCP3208_CS, LOW);

    SPI.transfer(comando);
    byte2 = SPI.transfer(byte2);
    byte3 = SPI.transfer(0x00);

    digitalWrite(MCP3208_CS, HIGH);

    uint16_t resultado;

    resultado = ((byte2 & 0x0F) << 8) | byte3;

    return resultado;
}

// =====================================================
// CONVERSÃO ADC -> VOLTS
// =====================================================

float adcParaVolts(uint16_t valor)
{
    return ((float)valor * MCP_VREF) / MCP_MAX;
}

// =====================================================
// AQUISIÇÃO SINCRONIZADA
// =====================================================

void adquirirAmostras()
{
    uint16_t anterior;
    uint16_t atual;

    // -----------------------------------------------
    // PROCURA CRUZAMENTO DE V1
    // -----------------------------------------------

    anterior = lerMCP3208(V1_CH);

    while (true)
    {
        atual = lerMCP3208(V1_CH);

        if (anterior < 2048 && atual >= 2048)
        {
            break;
        }

        anterior = atual;
    }

    // -----------------------------------------------
    // AQUISIÇÃO DOS 6 CANAIS
    // -----------------------------------------------

    for (int i = 0; i < SAMPLES; i++)
    {
        unsigned long t0 = micros();

        adc[0][i] = lerMCP3208(V1_CH);
        adc[1][i] = lerMCP3208(V2_CH);
        adc[2][i] = lerMCP3208(V3_CH);

        adc[3][i] = lerMCP3208(I1_CH);
        adc[4][i] = lerMCP3208(I2_CH);
        adc[5][i] = lerMCP3208(I3_CH);

        // Mantém o período de amostragem
        while ((micros() - t0) < Ts_us)
        {
        }
    }
}

// =====================================================
// REMOVE OFFSET
// =====================================================

void removerOffset(float offset[6])
{
    for (int c = 0; c < 6; c++)
    {
        double soma = 0;

        for (int i = 0; i < SAMPLES; i++)
        {
            soma += adc[c][i];
        }

        offset[c] = soma / SAMPLES;
    }
}

// =====================================================
// RMS + POTÊNCIA
// =====================================================

void calcularPotencias()
{
    float offset[6];

    removerOffset(offset);

    Ptotal = 0;
    Stotal = 0;

    for (int fase = 0; fase < 3; fase++)
    {
        double somaV2 = 0;
        double somaI2 = 0;
        double somaP = 0;

        for (int i = 0; i < SAMPLES; i++)
        {
            float v =
                ((float)adc[fase][i] - offset[fase]);

            float iCorr =
                ((float)adc[fase + 3][i] - offset[fase + 3]);

            // ADC -> tensão
            v *= MCP_VREF / MCP_MAX;
            iCorr *= MCP_VREF / MCP_MAX;

            somaV2 += v * v;
            somaI2 += iCorr * iCorr;
            somaP += v * iCorr;
        }

        Vrms[fase] = sqrt(somaV2 / SAMPLES);
        Irms[fase] = sqrt(somaI2 / SAMPLES);

        // Correção original do código
        Vrms[fase] *= 1.02541;

        PotAtiva[fase] = somaP / SAMPLES;

        PotAparente[fase] =
            Vrms[fase] * Irms[fase];

        if (PotAparente[fase] > 0.0001)
        {
            FP[fase] =
                PotAtiva[fase] / PotAparente[fase];
        }
        else
        {
            FP[fase] = 0;
        }

        Ptotal += PotAtiva[fase];
        Stotal += PotAparente[fase];
    }

    if (Stotal > 0.0001)
    {
        FPtotal = Ptotal / Stotal;
    }
    else
    {
        FPtotal = 0;
    }
}

// =====================================================
// CÁLCULO DA FFT
// =====================================================

void calcularFFT()
{
    float offset[6];

    removerOffset(offset);

    // -----------------------------------------------
    // FFT DO V1
    // -----------------------------------------------

    for (int i = 0; i < SAMPLES; i++)
    {
        vReal[i] =
            ((float)adc[0][i] - offset[0]);

        vImag[i] = 0;
    }

    // Janela Hamming
    FFT.windowing(
        vReal,
        SAMPLES,
        FFT_WIN_TYP_HAMMING,
        FFT_FORWARD
    );

    // FFT
    FFT.compute(
        vReal,
        vImag,
        SAMPLES,
        FFT_FORWARD
    );

    // -----------------------------------------------
    // ENCONTRA O MAIOR PICO
    // -----------------------------------------------

    int binMaior = 1;
    float maiorMag = 0;

    for (int i = 1; i < SAMPLES / 2; i++)
    {
        float mag =
            sqrt(
                vReal[i] * vReal[i] +
                vImag[i] * vImag[i]
            );

        if (mag > maiorMag)
        {
            maiorMag = mag;
            binMaior = i;
        }
    }

    binFundamental = binMaior;

    freqFundamental =
        ((float)binFundamental *
         SAMPLING_FREQUENCY) /
        SAMPLES;

    // -----------------------------------------------
    // FASE DOS 6 CANAIS
    // -----------------------------------------------

    for (int canal = 0; canal < 6; canal++)
    {
        for (int i = 0; i < SAMPLES; i++)
        {
            vReal[i] =
                ((float)adc[canal][i] -
                 offset[canal]);

            vImag[i] = 0;
        }

        FFT.windowing(
            vReal,
            SAMPLES,
            FFT_WIN_TYP_HAMMING,
            FFT_FORWARD
        );

        FFT.compute(
            vReal,
            vImag,
            SAMPLES,
            FFT_FORWARD
        );

        faseFund[canal] =
            atan2(
                vImag[binFundamental],
                vReal[binFundamental]
            )
            * 180.0 / PI;
    }

    // -----------------------------------------------
    // DEFASAGEM V x I
    // -----------------------------------------------

    phi[0] =
        faseFund[0] - faseFund[3];

    phi[1] =
        faseFund[1] - faseFund[4];

    phi[2] =
        faseFund[2] - faseFund[5];

    for (int f = 0; f < 3; f++)
    {
        while (phi[f] > 180)
            phi[f] -= 360;

        while (phi[f] < -180)
            phi[f] += 360;

        DPF[f] =
            cos(phi[f] * PI / 180.0);
    }
}

// =====================================================
// CÁLCULO DE HARMÔNICAS E THD
// =====================================================

void calcularTHD()
{
    float offset[6];

    removerOffset(offset);

    for (int fase = 0; fase < 3; fase++)
    {
        // =================================================
        // TENSÃO
        // =================================================

        for (int i = 0; i < SAMPLES; i++)
        {
            vReal[i] =
                ((float)adc[fase][i] -
                 offset[fase]);

            vImag[i] = 0;
        }

        FFT.windowing(
            vReal,
            SAMPLES,
            FFT_WIN_TYP_HAMMING,
            FFT_FORWARD
        );

        FFT.compute(
            vReal,
            vImag,
            SAMPLES,
            FFT_FORWARD
        );

        float somaH2_V = 0;
        float rmsTotal_V = 0;

        // H1 até H9
        for (int h = 1; h <= 9; h++)
        {
            int binH =
                binFundamental * h;

            if (binH >= SAMPLES / 2)
                break;

            float mag =
                sqrt(
                    vReal[binH] * vReal[binH] +
                    vImag[binH] * vImag[binH]
                );

            // Amplitude de pico
            float amp =
                (2.0 * mag) / SAMPLES;

            // Correção Hamming original
            amp *= 1.85;

            // ADC -> V
            amp *= MCP_VREF / MCP_MAX;

            // Correção original
            amp *= 1.02541;

            // Pico -> RMS
            float rms =
                amp / sqrt(2.0);

            harmBinV[fase][h] = binH;

            harmFreqV[fase][h] =
                (float)binH *
                SAMPLING_FREQUENCY /
                SAMPLES;

            harmAmpV[fase][h] = rms;

            if (h == 1)
            {
                VrmsFund[fase] = rms;
            }
            else
            {
                somaH2_V += rms * rms;
            }

            rmsTotal_V += rms * rms;
        }

        VrmsTotalFFT[fase] =
            sqrt(rmsTotal_V);

        if (VrmsFund[fase] > 0.005)
        {
            THDv[fase] =
                (
                    sqrt(somaH2_V) /
                    VrmsFund[fase]
                ) * 100.0;
        }
        else
        {
            THDv[fase] = 0;
        }

        // =================================================
        // CORRENTE
        // =================================================

        for (int i = 0; i < SAMPLES; i++)
        {
            vReal[i] =
                ((float)adc[fase + 3][i] -
                 offset[fase + 3]);

            vImag[i] = 0;
        }

        FFT.windowing(
            vReal,
            SAMPLES,
            FFT_WIN_TYP_HAMMING,
            FFT_FORWARD
        );

        FFT.compute(
            vReal,
            vImag,
            SAMPLES,
            FFT_FORWARD
        );

        float somaH2_I = 0;
        float rmsTotal_I = 0;

        for (int h = 1; h <= 9; h++)
        {
            int binH =
                binFundamental * h;

            if (binH >= SAMPLES / 2)
                break;

            float mag =
                sqrt(
                    vReal[binH] * vReal[binH] +
                    vImag[binH] * vImag[binH]
                );

            float amp =
                (2.0 * mag) / SAMPLES;

            // Correção Hamming
            amp *= 1.85;

            // ADC -> V
            amp *= MCP_VREF / MCP_MAX;

            // Pico -> RMS
            float rms =
                amp / sqrt(2.0);

            harmBinI[fase][h] = binH;

            harmFreqI[fase][h] =
                (float)binH *
                SAMPLING_FREQUENCY /
                SAMPLES;

            harmAmpI[fase][h] = rms;

            if (h == 1)
            {
                IrmsFund[fase] = rms;
            }
            else
            {
                somaH2_I += rms * rms;
            }

            rmsTotal_I += rms * rms;
        }

        IrmsTotalFFT[fase] =
            sqrt(rmsTotal_I);

        if (IrmsFund[fase] > 0.005)
        {
            THDi[fase] =
                (
                    sqrt(somaH2_I) /
                    IrmsFund[fase]
                ) * 100.0;
        }
        else
        {
            THDi[fase] = 0;
        }

        // =================================================
        // MAIOR HARMÔNICA DE TENSÃO
        // =================================================

        harmPicoV[fase] = 2;

        harmPicoAmpV[fase] =
            harmAmpV[fase][2];

        harmPicoFreqV[fase] =
            harmFreqV[fase][2];

        for (int h = 3; h <= 9; h++)
        {
            if (
                harmAmpV[fase][h] >
                harmPicoAmpV[fase]
            )
            {
                harmPicoAmpV[fase] =
                    harmAmpV[fase][h];

                harmPicoV[fase] = h;

                harmPicoFreqV[fase] =
                    harmFreqV[fase][h];
            }
        }

        // =================================================
        // MAIOR HARMÔNICA DE CORRENTE
        // =================================================

        harmPicoI[fase] = 2;

        harmPicoAmpI[fase] =
            harmAmpI[fase][2];

        harmPicoFreqI[fase] =
            harmFreqI[fase][2];

        for (int h = 3; h <= 9; h++)
        {
            if (
                harmAmpI[fase][h] >
                harmPicoAmpI[fase]
            )
            {
                harmPicoAmpI[fase] =
                    harmAmpI[fase][h];

                harmPicoI[fase] = h;

                harmPicoFreqI[fase] =
                    harmFreqI[fase][h];
            }
        }

        // =================================================
        // DF
        // =================================================

        if (fabs(DPF[fase]) > 0.001)
        {
            DF[fase] =
                FP[fase] / DPF[fase];
        }
        else
        {
            DF[fase] = 0;
        }
    }
}

// =====================================================
// IMPRESSÃO DE POTÊNCIAS
// =====================================================

void imprimirPotencias()
{
    Serial.println();
    Serial.println("========== POTENCIA ==========");

    for (int f = 0; f < 3; f++)
    {
        Serial.print("FASE ");
        Serial.println(f + 1);

        Serial.print("Vrms = ");
        Serial.println(Vrms[f], 4);

        Serial.print("Irms = ");
        Serial.println(Irms[f], 4);

        Serial.print("P = ");
        Serial.println(PotAtiva[f], 4);

        Serial.print("S = ");
        Serial.println(PotAparente[f], 4);

        Serial.print("FP = ");
        Serial.println(FP[f], 4);

        Serial.println();
    }

    Serial.println("----- TOTAL -----");

    Serial.print("Ptotal = ");
    Serial.println(Ptotal, 4);

    Serial.print("Stotal = ");
    Serial.println(Stotal, 4);

    Serial.print("FPtotal = ");
    Serial.println(FPtotal, 4);
}

// =====================================================
// IMPRESSÃO FFT
// =====================================================

void imprimirFFT()
{
    Serial.println();
    Serial.println("========== FFT ==========");

    Serial.print("Fundamental = ");
    Serial.print(freqFundamental, 2);
    Serial.println(" Hz");

    Serial.print("Bin Fundamental = ");
    Serial.println(binFundamental);

    Serial.print("Resolucao FFT = ");
    Serial.print(
        (float)SAMPLING_FREQUENCY /
        SAMPLES,
        2
    );
    Serial.println(" Hz/bin");

    Serial.println();

    for (int f = 0; f < 3; f++)
    {
        Serial.print("FASE ");
        Serial.println(f + 1);

        Serial.print("Defasagem = ");
        Serial.print(phi[f], 2);
        Serial.println(" graus");

        Serial.print("DPF = ");
        Serial.println(DPF[f], 4);

        Serial.println();
    }
}

// =====================================================
// IMPRESSÃO THD
// =====================================================

void imprimirTHD()
{
    Serial.println();
    Serial.println("========== THD E HARMONICAS ==========");

    for (int f = 0; f < 3; f++)
    {
        Serial.print("FASE ");
        Serial.println(f + 1);

        // ---------------------------------------------
        // TENSÃO
        // ---------------------------------------------

        Serial.println("--- TENSAO ---");

        Serial.print("Vrms Fundamental = ");
        Serial.println(
            VrmsFund[f] * GainV,
            4
        );

        Serial.print("Vrms Total FFT = ");
        Serial.println(
            VrmsTotalFFT[f] * GainV,
            4
        );

        Serial.print("THDv = ");
        Serial.print(THDv[f], 2);
        Serial.println(" %");

        Serial.println(
            "Espectro de Tensao (Vrms):"
        );

        for (int h = 1; h <= 9; h++)
        {
            Serial.print("  H");
            Serial.print(h);

            Serial.print(" (");
            Serial.print(
                harmFreqV[f][h],
                1
            );
            Serial.print(" Hz): ");

            Serial.println(
                harmAmpV[f][h] * GainV,
                4
            );
        }

        Serial.println();

        // ---------------------------------------------
        // CORRENTE
        // ---------------------------------------------

        Serial.println("--- CORRENTE ---");

        Serial.print("Irms Fundamental = ");
        Serial.println(
            IrmsFund[f] * GainI,
            4
        );

        Serial.print("Irms Total FFT = ");
        Serial.println(
            IrmsTotalFFT[f] * GainI,
            4
        );

        Serial.print("THDi = ");
        Serial.print(THDi[f], 2);
        Serial.println(" %");

        Serial.println(
            "Espectro de Corrente (Irms):"
        );

        for (int h = 1; h <= 9; h++)
        {
            Serial.print("  H");
            Serial.print(h);

            Serial.print(" (");
            Serial.print(
                harmFreqI[f][h],
                1
            );
            Serial.print(" Hz): ");

            Serial.println(
                harmAmpI[f][h] * GainI,
                4
            );
        }

        Serial.println();

        // ---------------------------------------------
        // FATOR DE POTÊNCIA
        // ---------------------------------------------

        Serial.print("DPF = ");
        Serial.println(DPF[f], 4);

        Serial.print("FP = ");
        Serial.println(FP[f], 4);

        Serial.print("DF = ");
        Serial.println(DF[f], 4);

        // ---------------------------------------------
        // MAIOR HARMÔNICA DE TENSÃO
        // ---------------------------------------------

        float porcentagemFundV = 0;

        if ((VrmsFund[f] * GainV) > 0.0001)
        {
            porcentagemFundV =
                (
                    (harmPicoAmpV[f] * GainV) /
                    (VrmsFund[f] * GainV)
                ) * 100.0;
        }

        Serial.print(
            "Maior Harmonica Tensao H"
        );

        Serial.print(
            harmPicoV[f]
        );

        Serial.print(" (");

        Serial.print(
            harmPicoFreqV[f],
            1
        );

        Serial.print(" Hz) Amplitude ");

        Serial.print(
            harmPicoAmpV[f] * GainV,
            4
        );

        Serial.print(" Vrms (");

        Serial.print(
            porcentagemFundV,
            1
        );

        Serial.println(
            "%) da fundamental"
        );

        // ---------------------------------------------
        // MAIOR HARMÔNICA DE CORRENTE
        // ---------------------------------------------

        float porcentagemFundI = 0;

        if ((IrmsFund[f] * GainI) > 0.0001)
        {
            porcentagemFundI =
                (
                    (harmPicoAmpI[f] * GainI) /
                    (IrmsFund[f] * GainI)
                ) * 100.0;
        }

        Serial.print(
            "Maior Harmonica Corrente H"
        );

        Serial.print(
            harmPicoI[f]
        );

        Serial.print(" (");

        Serial.print(
            harmPicoFreqI[f],
            1
        );

        Serial.print(" Hz) Amplitude ");

        Serial.print(
            harmPicoAmpI[f] * GainI,
            4
        );

        Serial.print(" Irms (");

        Serial.print(
            porcentagemFundI,
            1
        );

        Serial.println(
            "%) da fundamental"
        );

        Serial.println(
            "------------------------------------"
        );
    }
}

// =====================================================
// SETUP
// =====================================================

void setup()
{
    Serial.begin(9600);

    // -----------------------------------------------
    // CS DO MCP3208
    // -----------------------------------------------

    pinMode(MCP3208_CS, OUTPUT);
    digitalWrite(MCP3208_CS, HIGH);

    // -----------------------------------------------
    // SPI
    // -----------------------------------------------

    SPI.begin();

    // -----------------------------------------------
    // TEMPO DE AMOSTRAGEM
    // -----------------------------------------------

    Ts_us =
        round(
            1000000.0 /
            SAMPLING_FREQUENCY
        );

    Serial.println();
    Serial.println(
        "ANALISADOR TRIFASICO"
    );

    Serial.println(
        "ADC: MCP3208 - 12 BITS"
    );

    Serial.println(
        "VREF: 5.0 V"
    );

    Serial.print(
        "Fs = "
    );

    Serial.print(
        SAMPLING_FREQUENCY
    );

    Serial.println(
        " Hz"
    );

    Serial.print(
        "Ts = "
    );

    Serial.print(
        Ts_us
    );

    Serial.println(
        " us"
    );

    Serial.print(
        "Resolucao FFT = "
    );

    Serial.print(
        (float)SAMPLING_FREQUENCY /
        SAMPLES,
        2
    );

    Serial.println(
        " Hz/bin"
    );
}

// =====================================================
// LOOP
// =====================================================

void loop()
{
    adquirirAmostras();

    calcularPotencias();

    calcularFFT();

    calcularTHD();

    imprimirPotencias();

    imprimirFFT();

    imprimirTHD();

    delay(1000);
}