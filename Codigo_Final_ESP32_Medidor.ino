#include <WiFi.h>
#include <SPI.h>
#include <arduinoFFT.h>
#include <ESPmDNS.h>
#include <WebServer.h>
#include <Preferences.h>
#include <math.h>

// =====================================================
// PINAGEM SPI ESP32 (MCP3208)
// =====================================================
#define MCP3208_CS 5
#define SPI_SCK    18
#define SPI_MISO   19
#define SPI_MOSI   23

// =====================================================
// CONFIGURAÇÕES DE AMOSTRAGEM E ADC
// =====================================================
#define SAMPLES 256
#define SAMPLING_FREQUENCY 3840.0
#define MCP_VREF 5.0
#define MCP_MAX 4095.0
#define MEDIA_TAMANHO 10

#define SPI_CLOCK 2000000UL
SPISettings mcp3208SPI(SPI_CLOCK, MSBFIRST, SPI_MODE0);

// Fatores Base de Hardware
float baseGainV = 755.5556; 
float baseGainI = 400.0000; 
const float GainAmpV = 2.0000;

// Fatores de Calibração Individuais (CH0 a CH5)
float calibFactor[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
float erroPercentual[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

// Mapeamento Intercalado (CH0=V1, CH1=I1, CH2=V2, CH3=I2, CH4=V3, CH5=I3)
const uint8_t CH_V[3] = {0, 2, 4}; 
const uint8_t CH_I[3] = {1, 3, 5}; 

bool modoVLL = false;         
uint8_t modoSistema = 2; // 0 = Monofásico (CH0, CH1), 1 = Bifásico (CH0..CH3), 2 = Trifásico (CH0..CH5)

// Buffers e Variáveis Elétricas
uint16_t adc[6][SAMPLES];
float vReal[SAMPLES];
float vImag[SAMPLES];

float VrmsTrueTime[6]; // True RMS por canal bruto no tempo
float Vrms[3], Irms[3], PotAtiva[3], PotReativa[3], PotAparente[3], FP[3];
float VrmsMed[3], IrmsMed[3], PMed[3], QMed[3], SMed[3], FPMed[3];
float PtotalMed = 0, QtotalMed = 0, StotalMed = 0, FPtotalMed = 0;

// Registro de Picos Máximos e Mínimos por Canal (CH0..CH5)
float picoMaxCh[6] = {0, 0, 0, 0, 0, 0};
float picoMinCh[6] = {9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 9999.0};

double energiaAcumulada_kWh = 0.0;
unsigned long tempoAnterior_ms = 0;

float freqFundamental = 0;
float faseFund[6], phi[3], DPF[3], THDv[3], THDi[3];
float VrmsFund[3], IrmsFund[3], VrmsTotalFFT[3], IrmsTotalFFT[3], DF[3];
float THDvTotal = 0.0, THDiTotal = 0.0;

float harmAmpV[3][10], harmAmpI[3][10];
float harmPctV[3][10], harmPctI[3][10];
float timeWaveform[6][SAMPLES];

// Média Móvel
float bufVrms[3][MEDIA_TAMANHO], bufIrms[3][MEDIA_TAMANHO];
float bufP[3][MEDIA_TAMANHO], bufQ[3][MEDIA_TAMANHO], bufS[3][MEDIA_TAMANHO], bufFP[3][MEDIA_TAMANHO];
uint8_t idxMedia = 0;
bool mediaCheia = false;

unsigned long Ts_us;
bool wifiConectado = false;
WebServer server(80);
Preferences pref;

uint16_t lerMCP3208(uint8_t canal) {
    canal &= 0x07;
    uint8_t cmd1 = 0x06 | (canal >> 2);
    uint8_t cmd2 = (canal & 0x03) << 6;

    SPI.beginTransaction(mcp3208SPI);
    digitalWrite(MCP3208_CS, LOW);
    SPI.transfer(cmd1);
    uint8_t msb = SPI.transfer(cmd2);
    uint8_t lsb = SPI.transfer(0x00);
    digitalWrite(MCP3208_CS, HIGH);
    SPI.endTransaction();

    return ((msb & 0x0F) << 8) | lsb;
}

bool isTensaoChannel(uint8_t ch) {
    return (ch == 0 || ch == 2 || ch == 4);
}

void adquirirAmostras() {
    uint16_t anterior = lerMCP3208(0);
    uint16_t atual;
    unsigned long timeout = millis();

    while ((millis() - timeout) < 50) {
        atual = lerMCP3208(0);
        if (anterior < 2048 && atual >= 2048) break;
        anterior = atual;
    }

    for (int i = 0; i < SAMPLES; i++) {
        unsigned long t0 = micros();
        for (int c = 0; c < 6; c++) adc[c][i] = 0;

        if (modoSistema == 0) {
            adc[0][i] = lerMCP3208(0);
            adc[1][i] = lerMCP3208(1);
        } else if (modoSistema == 1) {
            adc[0][i] = lerMCP3208(0);
            adc[1][i] = lerMCP3208(1);
            adc[2][i] = lerMCP3208(2);
            adc[3][i] = lerMCP3208(3);
        } else {
            for (int c = 0; c < 6; c++) {
                adc[c][i] = lerMCP3208(c);
            }
        }
        while ((micros() - t0) < Ts_us);
    }
}

void processarSinais() {
    float offset[6] = {0};
    for (int c = 0; c < 6; c++) {
        double soma = 0;
        for (int i = 0; i < SAMPLES; i++) soma += adc[c][i];
        offset[c] = (float)(soma / SAMPLES);
    }

    for (int c = 0; c < 6; c++) {
        for (int i = 0; i < SAMPLES; i++) {
            timeWaveform[c][i] = ((float)adc[c][i] - offset[c]) * (MCP_VREF / MCP_MAX) * 1000.0;
        }
    }

    // FFT no CH0 para Frequência Fundamental
    for (int i = 0; i < SAMPLES; i++) {
        vReal[i] = ((float)adc[0][i] - offset[0]);
        vImag[i] = 0;
    }
    ArduinoFFT<float> FFT = ArduinoFFT<float>(vReal, vImag, SAMPLES, SAMPLING_FREQUENCY);
    FFT.windowing(FFTWindow::Hamming, FFTDirection::Forward);
    FFT.compute(FFTDirection::Forward);

    int binFundamental = 1;
    float maiorMag = 0;
    for (int i = 1; i < SAMPLES / 2; i++) {
        float mag = sqrt(vReal[i] * vReal[i] + vImag[i] * vImag[i]);
        if (mag > maiorMag) { maiorMag = mag; binFundamental = i; }
    }
    freqFundamental = ((float)binFundamental * SAMPLING_FREQUENCY) / SAMPLES;

    for (int ch = 0; ch < 6; ch++) {
        if (modoSistema == 0 && (ch != 0 && ch != 1)) continue;
        if (modoSistema == 1 && (ch == 4 || ch == 5)) continue;

        bool isTensao = isTensaoChannel(ch);
        uint8_t fase = ch / 2;

        // 1. CÁLCULO BRUTO TRUE RMS VIA SOMA QUADRÁTICA NO TEMPO
        double somaQuadrados = 0;
        for (int i = 0; i < SAMPLES; i++) {
            float valVoltsADC = ((float)adc[ch][i] - offset[ch]) * (MCP_VREF / MCP_MAX);
            somaQuadrados += (valVoltsADC * valVoltsADC);
        }
        float vRmsAdc = sqrt(somaQuadrados / SAMPLES);

        // Aplica ganho base e calibração para obter True RMS Real do Sinal
        float vRmsTrue = isTensao ? (((vRmsAdc / GainAmpV)) * baseGainV * calibFactor[ch])
                                  : (vRmsAdc * baseGainI * calibFactor[ch]);

        VrmsTrueTime[ch] = vRmsTrue;

        // Atualização dos Picos Máximo e Mínimo por Canal
        if (vRmsTrue > picoMaxCh[ch]) picoMaxCh[ch] = vRmsTrue;
        if (vRmsTrue < picoMinCh[ch] && vRmsTrue > 0.05) picoMinCh[ch] = vRmsTrue;

        if (isTensao) VrmsTotalFFT[fase] = vRmsTrue;
        else IrmsTotalFFT[fase] = vRmsTrue;

        // 2. FFT DAS FREQUÊNCIAS E DECOMPOSIÇÃO HARMÔNICA
        for (int i = 0; i < SAMPLES; i++) {
            vReal[i] = ((float)adc[ch][i] - offset[0]);
            vImag[i] = 0;
        }
        ArduinoFFT<float> FFT_C = ArduinoFFT<float>(vReal, vImag, SAMPLES, SAMPLING_FREQUENCY);
        FFT_C.windowing(FFTWindow::Hamming, FFTDirection::Forward);
        FFT_C.compute(FFTDirection::Forward);

        faseFund[ch] = atan2(vImag[binFundamental], vReal[binFundamental]) * 180.0 / M_PI;

        float somaH2 = 0;
        for (int h = 1; h <= 9; h++) {
            int binH = binFundamental * h;
            if (binH >= SAMPLES / 2) break;

            float mag = sqrt(vReal[binH] * vReal[binH] + vImag[binH] * vImag[binH]);
            
            float amp = (2.0 * mag / SAMPLES);
            if (isTensao) {
                amp *= 1.893;
                amp *= (MCP_VREF / MCP_MAX);
            } else {
                amp *= 1.87;
                amp *= (MCP_VREF / MCP_MAX);
            }
            
            float rms = isTensao ? (((amp / GainAmpV) / sqrt(2.0)) * baseGainV * calibFactor[ch]) 
                                 : ((amp / sqrt(2.0)) * baseGainI * calibFactor[ch]);

            if (isTensao) harmAmpV[fase][h] = rms;
            else harmAmpI[fase][h] = rms;

            if (h == 1) {
                if (isTensao) VrmsFund[fase] = rms;
                else IrmsFund[fase] = rms;
            } else {
                somaH2 += rms * rms;
            }
        }

        float fundVal = isTensao ? VrmsFund[fase] : IrmsFund[fase];
        for (int h = 1; h <= 9; h++) {
            float val = isTensao ? harmAmpV[fase][h] : harmAmpI[fase][h];
            float pct = (fundVal > 0.001) ? ((val / fundVal) * 100.0) : 0.0;
            if (isTensao) harmPctV[fase][h] = pct;
            else harmPctI[fase][h] = pct;
        }

        if (isTensao) {
            THDv[fase] = (VrmsFund[fase] > 0.1) ? ((sqrt(somaH2) / VrmsFund[fase]) * 100.0) : 0;
        } else {
            THDi[fase] = (IrmsFund[fase] > 0.01) ? ((sqrt(somaH2) / IrmsFund[fase]) * 100.0) : 0;
        }
    }

    float angCorr = modoVLL ? 30.0 : 0.0;

    for (int f = 0; f < 3; f++) {
        if ((modoSistema == 0 && f != 0) || (modoSistema == 1 && f == 2)) {
            Vrms[f] = 0; Irms[f] = 0; PotAtiva[f] = 0; PotReativa[f] = 0; PotAparente[f] = 0; FP[f] = 0;
            continue;
        }

        uint8_t chV = CH_V[f];
        uint8_t chI = CH_I[f];

        phi[f] = (faseFund[chV] - faseFund[chI]) - angCorr;
        while (phi[f] > 180.0) phi[f] -= 360.0;
        while (phi[f] < -180.0) phi[f] += 360.0;

        DPF[f] = cos(phi[f] * M_PI / 180.0);
        DF[f] = 1.0 / (sqrt(1.0 + pow(THDv[f]/100.0, 2)) * sqrt(1.0 + pow(THDi[f]/100.0, 2)));

        Irms[f] = IrmsTotalFFT[f];
        Vrms[f] = VrmsTotalFFT[f];

        PotAparente[f] = Vrms[f] * Irms[f];
        FP[f] = DPF[f] * DF[f];
        if (FP[f] > 1.0) FP[f] = 1.0;
        if (FP[f] < -1.0) FP[f] = -1.0;
        
        PotAtiva[f] = PotAparente[f] * FP[f];
        PotReativa[f] = PotAparente[f] * sin(phi[f] * M_PI / 180.0);
    }

    // THD Total do Sistema
    uint8_t numFasesAtivas = (modoSistema == 0) ? 1 : ((modoSistema == 1) ? 2 : 3);
    float somaTHDv2 = 0.0, somaTHDi2 = 0.0;
    for (int f = 0; f < numFasesAtivas; f++) {
        somaTHDv2 += pow(THDv[f], 2);
        somaTHDi2 += pow(THDi[f], 2);
    }
    THDvTotal = sqrt(somaTHDv2 / numFasesAtivas);
    THDiTotal = sqrt(somaTHDi2 / numFasesAtivas);
}

void atualizarMediaMovel() {
    for (int f = 0; f < 3; f++) {
        bufVrms[f][idxMedia] = Vrms[f];
        bufIrms[f][idxMedia] = Irms[f];
        bufP[f][idxMedia]    = PotAtiva[f];
        bufQ[f][idxMedia]    = PotReativa[f];
        bufS[f][idxMedia]    = PotAparente[f];
        bufFP[f][idxMedia]   = FP[f];
    }
    idxMedia++;
    if (idxMedia >= MEDIA_TAMANHO) { idxMedia = 0; mediaCheia = true; }

    uint8_t limite = mediaCheia ? MEDIA_TAMANHO : idxMedia;
    if (limite == 0) limite = 1;

    PtotalMed = 0; QtotalMed = 0; StotalMed = 0;
    for (int f = 0; f < 3; f++) {
        if ((modoSistema == 0 && f != 0) || (modoSistema == 1 && f == 2)) continue;

        double sV=0, sI=0, sP=0, sQ=0, sS=0, sFP=0;
        for (int k = 0; k < limite; k++) {
            sV += bufVrms[f][k]; sI += bufIrms[f][k]; sP += bufP[f][k]; 
            sQ += bufQ[f][k]; sS += bufS[f][k]; sFP += bufFP[f][k];
        }
        VrmsMed[f] = sV / limite; IrmsMed[f] = sI / limite;
        PMed[f] = sP / limite; QMed[f] = sQ / limite; SMed[f] = sS / limite; FPMed[f] = sFP / limite;

        PtotalMed += PMed[f]; QtotalMed += QMed[f]; StotalMed += SMed[f];
    }
    FPtotalMed = (StotalMed > 0.0001) ? (PtotalMed / StotalMed) : 0;

    unsigned long agora = millis();
    if (tempoAnterior_ms > 0) {
        double dt = (agora - tempoAnterior_ms) / 1000.0;
        energiaAcumulada_kWh += (PtotalMed / 1000.0) * (dt / 3600.0);
    }
    tempoAnterior_ms = agora;
}

void calibrarCanal(uint8_t ch, float valorReal) {
    if (ch >= 6 || valorReal <= 0) return;
    
    bool isTensao = isTensaoChannel(ch);
    uint8_t fase = ch / 2;
    float medidoAtual = isTensao ? VrmsMed[fase] : IrmsMed[fase];

    if (medidoAtual > 0.001) {
        float fatorNovo = calibFactor[ch] * (valorReal / medidoAtual);
        erroPercentual[ch] = ((medidoAtual - valorReal) / valorReal) * 100.0;
        calibFactor[ch] = fatorNovo;

        pref.begin("sys-config", false);
        pref.putFloat(("cal_" + String(ch)).c_str(), calibFactor[ch]);
        pref.end();
    }
}

void resetarCalibracao() {
    pref.begin("sys-config", false);
    for (int c = 0; c < 6; c++) {
        calibFactor[c] = 1.0;
        erroPercentual[c] = 0.0;
        pref.putFloat(("cal_" + String(c)).c_str(), 1.0);
    }
    pref.end();
}

void resetarPicos() {
    for (int c = 0; c < 6; c++) {
        picoMaxCh[c] = 0.0;
        picoMinCh[c] = 9999.0;
    }
}

// =====================================================
// SERVIDOR WEB E ENDPOINTS JSON
// =====================================================
void tratarDadosJSON() {
    String json = "{";
    json += "\"modo_sistema\":" + String(modoSistema) + ",";
    json += "\"modo_vll\":" + String(modoVLL ? 1 : 0) + ",";
    json += "\"freq_fund\":" + String(freqFundamental, 2) + ",";
    
    // Tensões e Correntes Médias Calibradas por Fase
    json += "\"vrms_fase\":[" + String(VrmsMed[0], 2) + "," + String(VrmsMed[1], 2) + "," + String(VrmsMed[2], 2) + "],";
    json += "\"irms_fase\":[" + String(IrmsMed[0], 2) + "," + String(IrmsMed[1], 2) + "," + String(IrmsMed[2], 2) + "],";

    // Potências e FP por Fase + Totais
    json += "\"p\":[" + String(PMed[0]/1000.0, 3) + "," + String(PMed[1]/1000.0, 3) + "," + String(PMed[2]/1000.0, 3) + "],";
    json += "\"q\":[" + String(QMed[0]/1000.0, 3) + "," + String(QMed[1]/1000.0, 3) + "," + String(QMed[2]/1000.0, 3) + "],";
    json += "\"s\":[" + String(SMed[0]/1000.0, 3) + "," + String(SMed[1]/1000.0, 3) + "," + String(SMed[2]/1000.0, 3) + "],";
    json += "\"fp\":[" + String(FPMed[0], 3) + "," + String(FPMed[1], 3) + "," + String(FPMed[2], 3) + "],";
    json += "\"ptotal_kw\":" + String(PtotalMed / 1000.0, 3) + ",";
    json += "\"qtotal_kvar\":" + String(QtotalMed / 1000.0, 3) + ",";
    json += "\"stotal_kva\":" + String(StotalMed / 1000.0, 3) + ",";
    json += "\"fptotal\":" + String(FPtotalMed, 3) + ",";
    json += "\"kwh\":" + String(energiaAcumulada_kWh, 4) + ",";
    
    // Ângulos, Defasagens e THD Total
    json += "\"fase_fund\":[";
    for(int i=0; i<6; i++) {
        json += String(faseFund[i], 1);
        if(i < 5) json += ",";
    }
    json += "],";
    
    json += "\"phi\":[" + String(phi[0], 1) + "," + String(phi[1], 1) + "," + String(phi[2], 1) + "],";
    json += "\"thd_v_total\":" + String(THDvTotal, 2) + ",";
    json += "\"thd_i_total\":" + String(THDiTotal, 2) + ",";
    
    // Canais, Picos e Harmônicos
    json += "\"channels\":[";
    for (uint8_t ch = 0; ch < 6; ch++) {
        bool isTensao = isTensaoChannel(ch);
        uint8_t fase = ch / 2;

        json += "{";
        json += "\"channel\":" + String(ch) + ",";
        json += "\"vpico_max\":" + String(picoMaxCh[ch], 2) + ",";
        json += "\"vpico_min\":" + String(picoMinCh[ch] > 9000.0 ? 0.0 : picoMinCh[ch], 2) + ",";
        json += "\"vrms_medido\":" + String(VrmsTrueTime[ch], 2) + ",";
        json += "\"ventrada_vrms\":" + String(isTensao ? VrmsMed[fase] : IrmsMed[fase], 2) + ",";
        json += "\"calib_factor\":" + String(calibFactor[ch], 4) + ",";
        json += "\"erro_pct\":" + String(erroPercentual[ch], 2) + ",";
        json += "\"thd\":" + String(isTensao ? THDv[fase] : THDi[fase], 2) + ",";

        json += "\"wave\":[";
        for (int i = 0; i < SAMPLES; i++) {
            json += String(timeWaveform[ch][i], 1);
            if (i < SAMPLES - 1) json += ",";
        }
        json += "],";

        json += "\"harmonics\":[";
        for (int h = 1; h <= 7; h++) {
            float valPct = isTensao ? harmPctV[fase][h] : harmPctI[fase][h];
            json += String(valPct, 2);
            if (h < 7) json += ",";
        }
        json += "],";

        json += "\"harmonics_rms\":[";
        for (int h = 1; h <= 7; h++) {
            float valRmsFreq = isTensao ? harmAmpV[fase][h] : harmAmpI[fase][h];
            json += String(valRmsFreq, 2);
            if (h < 7) json += ",";
        }
        json += "]";

        json += "}";
        if (ch < 5) json += ",";
    }
    json += "]}";
    server.send(200, "application/json", json);
}

void tratarCalibrarWeb() {
    if (server.hasArg("ch") && server.hasArg("val")) {
        uint8_t ch = server.arg("ch").toInt();
        float val = server.arg("val").toFloat();
        calibrarCanal(ch, val);
        server.send(200, "text/plain", "OK");
    } else {
        server.send(400, "text/plain", "Parametros invalidos");
    }
}

void tratarResetCalibWeb() {
    resetarCalibracao();
    server.send(200, "text/plain", "OK");
}

void tratarResetPicosWeb() {
    resetarPicos();
    server.send(200, "text/plain", "OK");
}

void tratarToggleVLLWeb() {
    if (server.hasArg("vll")) {
        modoVLL = (server.arg("vll").toInt() == 1);
        pref.begin("sys-config", false);
        pref.putBool("modo_vll", modoVLL);
        pref.end();
        server.send(200, "text/plain", "OK");
    } else {
        server.send(400, "text/plain", "Parametro vll ausente");
    }
}

void tratarSetModoSysWeb() {
    if (server.hasArg("modo")) {
        uint8_t m = server.arg("modo").toInt();
        if (m <= 2) {
            modoSistema = m;
            pref.begin("sys-config", false);
            pref.putUChar("modo_sys", modoSistema);
            pref.end();
            server.send(200, "text/plain", "OK");
            return;
        }
    }
    server.send(400, "text/plain", "Modo invalido");
}

void tratarRaizWeb() {
    String html = R"rawliteral(
<!DOCTYPE html>
<html lang="pt-BR">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Analisador de Energia ESP32</title>
    <style>
        * { box-sizing: border-box; margin: 0; padding: 0; font-family: 'Segoe UI', sans-serif; }
        body { background-color: #0b0e14; color: #e1e6ed; padding: 15px; }
        
        /* GUIAS ESTILO NAVEGADOR */
        .tab-bar { display: flex; background: #161b22; border-bottom: 2px solid #30363d; border-radius: 8px 8px 0 0; padding: 5px 5px 0 5px; gap: 4px; overflow-x: auto; }
        .tab-nav { padding: 10px 20px; background: #21262d; color: #8b949e; border: 1px solid #30363d; border-bottom: none; border-radius: 6px 6px 0 0; cursor: pointer; font-weight: 600; font-size: 0.9rem; transition: 0.2s; white-space: nowrap; }
        .tab-nav:hover { background: #30363d; color: #f0f6fc; }
        .tab-nav.active { background: #0d1117; color: #58a6ff; border-top: 2px solid #58a6ff; font-size: 0.95rem; }
        
        .tab-content { display: none; background: #0d1117; border: 1px solid #30363d; border-top: none; padding: 20px; border-radius: 0 0 8px 8px; }
        .tab-content.active { display: block; }

        /* CARDS E CONTROLES */
        .controls-bar { display: flex; flex-wrap: wrap; gap: 15px; margin-bottom: 20px; background: #161b22; padding: 12px; border-radius: 8px; align-items: center; border: 1px solid #30363d; }
        .dashboard-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 15px; margin-bottom: 20px; }
        .card-metric { background: #161b22; border-radius: 8px; padding: 15px; text-align: center; border: 1px solid #30363d; }
        .metric-title { font-size: 0.75rem; text-transform: uppercase; font-weight: bold; margin-bottom: 5px; letter-spacing: 0.5px; }
        .metric-val { font-size: 1.3rem; font-weight: bold; }
        
        .t-green { color: #3fb950; } .t-cyan { color: #58a6ff; }
        .t-purple { color: #bc8cff; } .t-yellow { color: #d29922; }
        .t-orange { color: #f0883e; } .t-red { color: #f85149; }

        .card-calib { background: #161b22; border-radius: 8px; padding: 15px; border: 1px solid #f0883e; display: flex; flex-direction: column; align-items: center; justify-content: center; }
        .calib-input-group { display: flex; gap: 8px; margin-top: 8px; }
        .calib-input { background: #0d1117; border: 1px solid #30363d; color: #fff; padding: 6px; border-radius: 4px; width: 100px; text-align: center; }

        /* TABELAS */
        table { width: 100%; border-collapse: collapse; margin-top: 15px; background: #161b22; border-radius: 8px; overflow: hidden; }
        th, td { padding: 12px; text-align: center; border-bottom: 1px solid #30363d; }
        th { background: #21262d; color: #8b949e; text-transform: uppercase; font-size: 0.8rem; }

        canvas { width: 100%; height: 260px; background-color: #080a0f; border-radius: 6px; border: 1px solid #30363d; margin-top: 10px; }
        .btn-custom { background: #21262d; color: #58a6ff; border: 1px solid #30363d; padding: 8px 15px; border-radius: 6px; cursor: pointer; font-weight: bold; }
        .btn-custom:hover { background: #30363d; }
        .btn-custom.active { background: #58a6ff; color: #000; }
    </style>
</head>
<body>

    <!-- GUIAS NAVEGADOR -->
    <div class="tab-bar">
        <button class="tab-nav active" onclick="abrirGuia(event, 'gGeral')">🌐 Visão Geral</button>
        <button class="tab-nav" onclick="abrirGuia(event, 'gCanaisFP')">🎛️ Canais e Calibração</button>
        <button class="tab-nav" onclick="abrirGuia(event, 'gTenCorDefasagem')">⚡ Tensão, Corrente e Defasagem</button>
        <button class="tab-nav" onclick="abrirGuia(event, 'gPotencias')">📊 Potências P, Q, S, FP</button>
        <button class="tab-nav" onclick="abrirGuia(event, 'gHarmonicos')">📈 Análise Harmônica</button>
    </div>

    <!-- CONTROLES FIXOS -->
    <div class="controls-bar">
        <label>Modo do Sistema:</label>
        <select id="selectModoSys" onchange="alterarModoSistema(this.value)" style="background:#21262d; color:#fff; border:1px solid #30363d; padding:6px; border-radius:4px;">
            <option value="0">Monofásico (CH0, CH1)</option>
            <option value="1">Bifásico (CH0..CH3)</option>
            <option value="2">Trifásico (CH0..CH5)</option>
        </select>
        <button class="btn-custom" id="btnModeV" onclick="alternarModoVLL()">Vff / Vf</button>
        <button class="btn-custom" style="color:#f85149;" onclick="resetarCalibracaoWeb()">Resetar Calibração RMS</button>
        <button class="btn-custom" style="color:#f0883e;" onclick="resetarPicosWeb()">Resetar Máx/Mín</button>
        <span style="margin-left:auto; font-size:0.85rem; color:#8b949e;" id="freqTxt">Freq: 60.0 Hz</span>
    </div>

    <!-- GUIA 1: VISÃO GERAL -->
    <div id="gGeral" class="tab-content active">
        <div class="dashboard-grid">
            <div class="card-metric"><div class="metric-title t-green">Ativa Total (P)</div><div class="metric-val t-green" id="pTotG">0.00 kW</div></div>
            <div class="card-metric"><div class="metric-title t-purple">Reativa Total (Q)</div><div class="metric-val t-purple" id="qTotG">0.00 kVAr</div></div>
            <div class="card-metric"><div class="metric-title t-yellow">Aparente Total (S)</div><div class="metric-val t-yellow" id="sTotG">0.00 kVA</div></div>
            <div class="card-metric"><div class="metric-title t-cyan">FP Total</div><div class="metric-val t-cyan" id="fpTotG">0.000</div></div>
            <div class="card-metric"><div class="metric-title t-orange">%THD V / %THD I (Total)</div><div class="metric-val t-orange" id="thdG">0.0% / 0.0%</div></div>
            <div class="card-metric"><div class="metric-title t-green">Energia Acumulada</div><div class="metric-val t-green" id="kwhG">0.0000 kWh</div></div>
        </div>
        <h3 style="margin-top:15px;">Osciloscópio Geral (Todos os Sinais)</h3>
        <canvas id="canvasGeral"></canvas>
    </div>

    <!-- GUIA 2: CANAIS E CALIBRAÇÃO (COM MÍNIMO E RESET DE PICOS) -->
    <div id="gCanaisFP" class="tab-content">
        <div style="display:flex; gap:10px; margin-bottom:15px;" id="btnCanaisNav"></div>
        <div class="dashboard-grid">
            <div class="card-metric"><div class="metric-title t-cyan">Sinal Bruto RMS</div><div class="metric-val t-cyan" id="cRms">0.00</div></div>
            <div class="card-metric"><div class="metric-title t-yellow">Calibrado RMS</div><div class="metric-val t-yellow" id="cCalib">0.00</div></div>
            <div class="card-metric"><div class="metric-title t-orange">Fator / Erro %</div><div class="metric-val t-orange" id="cInfo">1.000 (0.0%)</div></div>
            <div class="card-metric"><div class="metric-title t-red">Pico Máximo</div><div class="metric-val t-red" id="cPicoMax">0.00</div></div>
            <div class="card-metric"><div class="metric-title t-green">Pico Mínimo</div><div class="metric-val t-green" id="cPicoMin">0.00</div></div>
            <div class="card-calib">
                <div class="metric-title t-orange" id="calibTitle">Calibrar Canal Atual</div>
                <div class="calib-input-group">
                    <input type="number" step="0.1" id="valMultimetro" class="calib-input" placeholder="Valor Real">
                    <button class="btn-custom" onclick="executarCalibracao()">Calibrar</button>
                </div>
            </div>
        </div>
        <h3>Onda do Canal Selecionado</h3>
        <canvas id="canvasCanal"></canvas>
    </div>

    <!-- GUIA 3: TENSÃO, CORRENTE E DEFASAGEM (COM VALORES DE V E I) -->
    <div id="gTenCorDefasagem" class="tab-content">
        <h3>Defasagem Angular e Casamento V x I</h3>
        <table>
            <thead>
                <tr>
                    <th>Par Comparado</th>
                    <th>Tensão RMS (V)</th>
                    <th>Corrente RMS (A)</th>
                    <th>Defasagem V x I (ϕ)</th>
                    <th>Ângulo Tensao V(f)</th>
                    <th>Ângulo Corrente I(f)</th>
                    <th>DPF (Fator Deslocamento)</th>
                </tr>
            </thead>
            <tbody id="tblDefasagem"></tbody>
        </table>
        <h3 style="margin-top:20px;">Formas de Onda Sobrepostas (Tensão x Corrente)</h3>
        <canvas id="canvasVX"></canvas>
    </div>

    <!-- GUIA 4: POTÊNCIAS POR FASE E TOTAL -->
    <div id="gPotencias" class="tab-content">
        <h3>Balanço de Potências e FP do Sistema</h3>
        <table>
            <thead>
                <tr>
                    <th>Fase / Grandeza</th>
                    <th>Potência Ativa (P)</th>
                    <th>Potência Reativa (Q)</th>
                    <th>Potência Aparente (S)</th>
                    <th>Fator de Potência (FP)</th>
                </tr>
            </thead>
            <tbody id="tblPotencias"></tbody>
        </table>
    </div>

    <!-- GUIA 5: ANÁLISE HARMÔNICA -->
    <div id="gHarmonicos" class="tab-content">
        <h3>Espectro Harmónico por Canal</h3>
        <div style="display:flex; gap:10px; margin-bottom:15px;" id="btnHarmNav"></div>
        <canvas id="canvasHarm"></canvas>
        <h3 style="margin-top:20px;">Tabela Detalhada das Harmônicas (1ª até 7ª)</h3>
        <table>
            <thead>
                <tr>
                    <th>Ordem Harmônica</th>
                    <th>Amplitude RMS</th>
                    <th>Frequência (Hz)</th>
                    <th>Proporção % da Fundamental</th>
                </tr>
            </thead>
            <tbody id="tblHarmonicos"></tbody>
        </table>
    </div>

    <script>
        let canalAtual = 0;
        let harmCanalAtual = 0;
        let dadosGlobais = null;
        let isVLL = false;
        const nomesCanais = ["CH0 (V1)", "CH1 (I1)", "CH2 (V2)", "CH3 (I2)", "CH4 (V3)", "CH5 (I3)"];
        const cores = ["#ff3d00", "#ff9100", "#00e676", "#00e5ff", "#29b6f6", "#e040fb"];

        function abrirGuia(evt, idGuia) {
            document.querySelectorAll('.tab-content').forEach(c => c.classList.remove('active'));
            document.querySelectorAll('.tab-nav').forEach(b => b.classList.remove('active'));
            document.getElementById(idGuia).classList.add('active');
            evt.currentTarget.classList.add('active');
            if(dadosGlobais) renderizarInterface();
        }

        function alterarModoSistema(modo) {
            fetch(`/set_modo_sys?modo=${modo}`).then(() => atualizarDados());
        }

        function alternarModoVLL() {
            fetch(`/toggle_vll?vll=${isVLL ? 0 : 1}`).then(() => atualizarDados());
        }

        function resetarCalibracaoWeb() {
            if (confirm("Resetar calibração de todos os canais para 1.0?")) {
                fetch('/reset_calib').then(() => atualizarDados());
            }
        }

        function resetarPicosWeb() {
            fetch('/reset_picos').then(() => atualizarDados());
        }

        function executarCalibracao() {
            const val = document.getElementById('valMultimetro').value;
            if (!val || val <= 0) return alert("Digite um valor válido!");
            fetch(`/calibrar?ch=${canalAtual}&val=${val}`).then(() => {
                document.getElementById('valMultimetro').value = '';
                atualizarDados();
            });
        }

        function renderizarInterface() {
            const p = dadosGlobais;
            isVLL = p.modo_vll === 1;
            document.getElementById('selectModoSys').value = p.modo_sistema;
            document.getElementById('btnModeV').innerText = isVLL ? 'Modo: Vll (+30°)' : 'Modo: Vf';
            document.getElementById('freqTxt').innerText = `Freq: ${p.freq_fund.toFixed(1)} Hz`;
            
            // GUIA 1: VISÃO GERAL
            document.getElementById('pTotG').innerText = p.ptotal_kw.toFixed(3) + ' kW';
            document.getElementById('qTotG').innerText = p.qtotal_kvar.toFixed(3) + ' kVAr';
            document.getElementById('sTotG').innerText = p.stotal_kva.toFixed(3) + ' kVA';
            document.getElementById('fpTotG').innerText = p.fptotal.toFixed(3);
            document.getElementById('thdG').innerText = `${p.thd_v_total.toFixed(1)}% / ${p.thd_i_total.toFixed(1)}%`;
            document.getElementById('kwhG').innerText = p.kwh.toFixed(4) + ' kWh';
            desenharGeral(p);

            // GUIA 2: CANAIS (COM MÁXIMO E MÍNIMO)
            let htmlBtns = '';
            const maxC = p.modo_sistema === 0 ? 2 : (p.modo_sistema === 1 ? 4 : 6);
            if (canalAtual >= maxC) canalAtual = 0;
            for(let i=0; i<maxC; i++) {
                htmlBtns += `<button class="btn-custom ${i===canalAtual?'active':''}" onclick="canalAtual=${i};renderizarInterface()">${nomesCanais[i]}</button>`;
            }
            document.getElementById('btnCanaisNav').innerHTML = htmlBtns;
            const c = p.channels[canalAtual];
            const un = (canalAtual % 2 === 0) ? 'V' : 'A';
            document.getElementById('cRms').innerText = c.vrms_medido.toFixed(2) + ' ' + un;
            document.getElementById('cCalib').innerText = c.ventrada_vrms.toFixed(2) + ' ' + un;
            document.getElementById('cInfo').innerText = `${c.calib_factor.toFixed(3)} (${c.erro_pct.toFixed(1)}%)`;
            document.getElementById('cPicoMax').innerText = c.vpico_max.toFixed(2) + ' ' + un;
            document.getElementById('cPicoMin').innerText = c.vpico_min.toFixed(2) + ' ' + un;
            document.getElementById('calibTitle').innerText = `Calibrar ${nomesCanais[canalAtual]}`;
            desenharOnda('canvasCanal', [c.wave], [cores[canalAtual]]);

            // GUIA 3: DEFASAGEM E VALORES RMS DE TENSÃO E CORRENTE
            let htmlDef = '';
            const fases = p.modo_sistema === 0 ? 1 : (p.modo_sistema === 1 ? 2 : 3);
            for(let f=0; f<fases; f++) {
                const dpf = Math.cos(p.phi[f] * Math.PI / 180.0);
                htmlDef += `<tr>
                    <td>Fase ${String.fromCharCode(65+f)} (CH${f*2} x CH${f*2+1})</td>
                    <td class="t-yellow">${p.vrms_fase[f].toFixed(2)} V</td>
                    <td class="t-cyan">${p.irms_fase[f].toFixed(2)} A</td>
                    <td class="t-orange">${p.phi[f].toFixed(1)}°</td>
                    <td>${p.fase_fund[f*2].toFixed(1)}°</td>
                    <td>${p.fase_fund[f*2+1].toFixed(1)}°</td>
                    <td class="t-cyan">${dpf.toFixed(3)}</td>
                </tr>`;
            }
            document.getElementById('tblDefasagem').innerHTML = htmlDef;
            desenharSobrepostaVX(p);

            // GUIA 4: POTÊNCIAS
            let htmlPot = '';
            for(let f=0; f<fases; f++) {
                htmlPot += `<tr>
                    <td>Fase ${String.fromCharCode(65+f)}</td>
                    <td class="t-green">${p.p[f].toFixed(3)} kW</td>
                    <td class="t-purple">${p.q[f].toFixed(3)} kVAr</td>
                    <td class="t-yellow">${p.s[f].toFixed(3)} kVA</td>
                    <td class="t-cyan">${p.fp[f].toFixed(3)}</td>
                </tr>`;
            }
            htmlPot += `<tr style="font-weight:bold; background:#21262d;">
                <td>TOTAL SISTEMA</td>
                <td class="t-green">${p.ptotal_kw.toFixed(3)} kW</td>
                <td class="t-purple">${p.qtotal_kvar.toFixed(3)} kVAr</td>
                <td class="t-yellow">${p.stotal_kva.toFixed(3)} kVA</td>
                <td class="t-cyan">${p.fptotal.toFixed(3)}</td>
            </tr>`;
            document.getElementById('tblPotencias').innerHTML = htmlPot;

            // GUIA 5: HARMÔNICOS
            let htmlHBtns = '';
            if (harmCanalAtual >= maxC) harmCanalAtual = 0;
            for(let i=0; i<maxC; i++) {
                htmlHBtns += `<button class="btn-custom ${i===harmCanalAtual?'active':''}" onclick="harmCanalAtual=${i};renderizarInterface()">${nomesCanais[i]}</button>`;
            }
            document.getElementById('btnHarmNav').innerHTML = htmlHBtns;
            
            const hc = p.channels[harmCanalAtual];
            const isV = harmCanalAtual % 2 === 0;
            let htmlH = '';
            for(let h=1; h<=7; h++) {
                const freqH = p.freq_fund * h;
                const ampRms = hc.harmonics_rms[h-1];
                const pct = hc.harmonics[h-1];
                htmlH += `<tr>
                    <td>${h===1 ? '1ª (Fundamental)' : h + 'ª Harmônica'}</td>
                    <td class="t-yellow">${ampRms.toFixed(2)} ${isV?'V':'A'}</td>
                    <td>${freqH.toFixed(1)} Hz</td>
                    <td class="t-cyan">${pct.toFixed(2)} %</td>
                </tr>`;
            }
            document.getElementById('tblHarmonicos').innerHTML = htmlH;
            desenharFFTBarras(hc.harmonics);
        }

        function desenharOnda(canvasId, arraysOndas, arrayCores) {
            const canvas = document.getElementById(canvasId);
            const ctx = canvas.getContext('2d');
            const w = canvas.width = canvas.clientWidth, h = canvas.height = canvas.clientHeight;
            ctx.clearRect(0, 0, w, h);
            
            ctx.strokeStyle = '#30363d'; ctx.lineWidth = 1; ctx.beginPath();
            ctx.moveTo(0, h/2); ctx.lineTo(w, h/2); ctx.stroke();

            arraysOndas.forEach((wave, idx) => {
                if(!wave || wave.length === 0) return;
                let max = Math.max(...wave.map(Math.abs)); if(max<1) max=1;
                ctx.strokeStyle = arrayCores[idx]; ctx.lineWidth = 2; ctx.beginPath();
                const step = w / (wave.length - 1);
                for(let i=0; i<wave.length; i++) {
                    const x = i * step, y = (h/2) - (wave[i]/max)*(h/2 - 15);
                    if(i===0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
                }
                ctx.stroke();
            });
        }

        function desenharGeral(p) {
            const maxC = p.modo_sistema === 0 ? 2 : (p.modo_sistema === 1 ? 4 : 6);
            let waves = [], cList = [];
            for(let i=0; i<maxC; i++) { waves.push(p.channels[i].wave); cList.push(cores[i]); }
            desenharOnda('canvasGeral', waves, cList);
        }

        function desenharSobrepostaVX(p) {
            const maxC = p.modo_sistema === 0 ? 2 : (p.modo_sistema === 1 ? 4 : 6);
            let waves = [], cList = [];
            for(let i=0; i<maxC; i++) { waves.push(p.channels[i].wave); cList.push(cores[i]); }
            desenharOnda('canvasVX', waves, cList);
        }

        function desenharFFTBarras(harmonics) {
            const canvas = document.getElementById('canvasHarm');
            const ctx = canvas.getContext('2d');
            const w = canvas.width = canvas.clientWidth, h = canvas.height = canvas.clientHeight;
            ctx.clearRect(0, 0, w, h);
            const gap = 20, num = harmonics.length, barW = (w - gap * (num + 1)) / num;

            for(let i=0; i<num; i++) {
                const pct = harmonics[i];
                const barH = (pct / 100.0) * (h - 45);
                const x = gap + i * (barW + gap), y = h - barH - 25;
                ctx.fillStyle = '#58a6ff'; ctx.fillRect(x, y, barW, barH);
                ctx.fillStyle = '#8b949e'; ctx.font = '12px sans-serif'; ctx.textAlign = 'center';
                ctx.fillText(`${i+1}ª`, x + barW/2, h - 8);
                ctx.fillText(`${pct.toFixed(1)}%`, x + barW/2, y - 5);
            }
        }

        function atualizarDados() {
            fetch('/data').then(r => r.json()).then(payload => {
                dadosGlobais = payload;
                renderizarInterface();
            });
        }

        setInterval(atualizarDados, 1000);
        atualizarDados();
    </script>
</body>
</html>
)rawliteral";
    server.send(200, "text/html", html);
}

void escanearRedesWiFi() {
    Serial.println("SCAN_START");
    WiFi.mode(WIFI_STA);
    WiFi.disconnect();
    delay(100);
    int n = WiFi.scanNetworks();
    if (n == 0) Serial.println("NO_NETWORKS");
    else {
        for (int i = 0; i < n; ++i) {
            Serial.print("WIFI:"); Serial.print(WiFi.SSID(i)); Serial.print(";");
            Serial.print(WiFi.RSSI(i)); Serial.print(";");
            Serial.println(WiFi.encryptionType(i) == WIFI_AUTH_OPEN ? "ABERTA" : "PROTEGIDA");
            delay(10);
        }
    }
    Serial.println("SCAN_END");
}

bool conectarWiFi(String ssid, String senha, String nome) {
    nome.replace("_", "-"); nome.replace(" ", ""); nome.toLowerCase();
    WiFi.disconnect(true, true); delay(200);
    WiFi.setHostname(nome.c_str());
    WiFi.mode(WIFI_STA);
    WiFi.setSleep(false);
    WiFi.begin(ssid.c_str(), senha.c_str());

    int tentativas = 0;
    while (WiFi.status() != WL_CONNECTED && tentativas < 25) { delay(500); tentativas++; }

    if (WiFi.status() == WL_CONNECTED) {
        wifiConectado = true;
        if (MDNS.begin(nome.c_str())) MDNS.addService("http", "tcp", 80);
        server.on("/", tratarRaizWeb);
        server.on("/data", tratarDadosJSON);
        server.on("/calibrar", tratarCalibrarWeb);
        server.on("/reset_calib", tratarResetCalibWeb);
        server.on("/reset_picos", tratarResetPicosWeb);
        server.on("/toggle_vll", tratarToggleVLLWeb);
        server.on("/set_modo_sys", tratarSetModoSysWeb);
        server.begin();
        Serial.print("CONNECTED:"); Serial.println(WiFi.localIP());
        return true;
    }
    Serial.println("CONNECT_FAILED");
    return false;
}

void setup() {
    Serial.begin(115200);
    pinMode(MCP3208_CS, OUTPUT);
    digitalWrite(MCP3208_CS, HIGH);

    SPI.begin(SPI_SCK, SPI_MISO, SPI_MOSI, MCP3208_CS);
    Ts_us = round(1000000.0 / SAMPLING_FREQUENCY);

    pref.begin("sys-config", false);
    modoSistema = pref.getUChar("modo_sys", 2);
    modoVLL = pref.getBool("modo_vll", false);
    for (int c = 0; c < 6; c++) {
        calibFactor[c] = pref.getFloat(("cal_" + String(c)).c_str(), 1.0);
    }
    String savedSSID = pref.getString("ssid", "");
    String savedSenha = pref.getString("senha", "");
    String savedNome = pref.getString("nome", "j42");
    pref.end();

    if (savedSSID.length() > 0) conectarWiFi(savedSSID, savedSenha, savedNome);
}

void loop() {
    if (Serial.available() > 0) {
        String cmd = Serial.readStringUntil('\n'); cmd.trim();
        if (cmd == "SCAN") {
            escanearRedesWiFi();
        } else if (cmd.startsWith("SET_CONFIG:")) {
            String param = cmd.substring(11);
            int idx = param.indexOf(';');
            if (idx != -1) {
                modoSistema = param.substring(0, idx).toInt();
                modoVLL = (param.substring(idx + 1).toInt() == 1);

                pref.begin("sys-config", false);
                pref.putUChar("modo_sys", modoSistema);
                pref.putBool("modo_vll", modoVLL);
                pref.end();
                Serial.println("CONFIG_OK");
            }
        } else if (cmd.startsWith("CONNECT:")) {
            String dados = cmd.substring(8);
            int idx1 = dados.indexOf(';');
            int idx2 = dados.indexOf(';', idx1 + 1);
            if (idx1 != -1) {
                String ssid = dados.substring(0, idx1);
                String senha = (idx2 != -1) ? dados.substring(idx1 + 1, idx2) : dados.substring(idx1 + 1);
                String nome = (idx2 != -1) ? dados.substring(idx2 + 1) : "j42";

                if (conectarWiFi(ssid, senha, nome)) {
                    pref.begin("sys-config", false);
                    pref.putString("ssid", ssid);
                    pref.putString("senha", senha);
                    pref.putString("nome", nome);
                    pref.end();
                }
            }
        }
    }

    adquirirAmostras();
    processarSinais();
    atualizarMediaMovel();

    if (wifiConectado) server.handleClient();
}