/**
 * Visualisation of power flow
 *
 * Copyright (c) 2011 Steven Blair
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */


final int MAX_WIDTH = 1100;
final int MAX_HEIGHT = 650;
final int X_START = 350;
final int Y_START = 50;
final int IMPEDANCE_START_Y = 75;
final int PHASORS_START_X = 300;
final int WAVEFORMS_START_X = PHASORS_START_X + 100;
final int WAVEFORMS_END_X = 800;
final int ONE_PU_HEIGHT = 40;

final int Y_INTERVAL = 100;
final int VA_START = 220;
final int VB_START = VA_START + Y_INTERVAL;
final int I_START = VB_START + Y_INTERVAL;
final int P_START = I_START + Y_INTERVAL;
final int Q_START = P_START + Y_INTERVAL;

final float MOUSE_OVER_LINE_DISTANCE_THRESHOLD = 5.0;
final int OFFSET_ORIGINAL = 3;

final color VaColor = color(180, 33, 38);	// red
final color VbColor = color(226, 99, 102);	// salmon
final color IColor = color(34, 177, 76);	// green
final color PColor = color(36, 78, 198);	// blue
final color QColor = color(200, 191, 231);	// blueish

// initial values defined in JavaScript code
float VaMag;
float VaPhase;
float VbMag;
float VbPhase;
float R;
float X;
float f;

// derived values
float w;        // system angular speed
float C;
float L;
float Zmag;
float meanP, meanQ, rmsI, rmsdV;

float magA = VaMag;
float magB = VbMag;
float magC = R;
float phaseA = VaPhase;
float phaseB = VbPhase;
float phaseC = X;

float Ts = 0.00005;
float totalTime = 0.060;
int ITERATIONS = int(totalTime / Ts);

float Va[] = new float[ITERATIONS];
float Vb[] = new float[ITERATIONS];
float dV[] = new float[ITERATIONS];
float i[] = new float[ITERATIONS];
float P[] = new float[ITERATIONS];
float Q[] = new float[ITERATIONS];


boolean mouseIsOverLine(float x1, float y1, float x2, float y2) {
  float d = dist(x1, y1, x2, y2);
  float d1 = dist(x1, y1, mouseX, mouseY);
  float d2 = dist(x2, y2, mouseX, mouseY);

  // distance between vertices must be similar to sum of distances from each vertex to mouse
  if (d1 + d2 < d + MOUSE_OVER_LINE_DISTANCE_THRESHOLD) {
    return true;
  }

  return false;
}

void drawPowerSystem() {
  int centreX = WAVEFORMS_START_X +(WAVEFORMS_END_X - WAVEFORMS_START_X) / 2;
  int lineLengths = 130;
  int terminalHalfHeight = 15;
  int impedanceWidth = 60;
  int currentArraySize = 7;
  fill(0);
  stroke(255);

  // lines out from impedance
  line(centreX - lineLengths, IMPEDANCE_START_Y, centreX + lineLengths, IMPEDANCE_START_Y);

  // terminals
  line(centreX - lineLengths, IMPEDANCE_START_Y - terminalHalfHeight, centreX - lineLengths, IMPEDANCE_START_Y + terminalHalfHeight);
  line(centreX + lineLengths, IMPEDANCE_START_Y - terminalHalfHeight, centreX + lineLengths, IMPEDANCE_START_Y + terminalHalfHeight);

  // impedance
  rectMode(CENTER);
  rect(centreX, IMPEDANCE_START_Y, impedanceWidth, 18);

  // current arrow
  stroke(IColor);
  strokeWeight(3)
    line((centreX - lineLengths) + ((centreX - impedanceWidth / 2)  - (centreX - lineLengths)) / 2, IMPEDANCE_START_Y, (centreX - lineLengths) + ((centreX - impedanceWidth / 2)  - (centreX - lineLengths)) / 2 - currentArraySize, IMPEDANCE_START_Y + currentArraySize);
  line((centreX - lineLengths) + ((centreX - impedanceWidth / 2)  - (centreX - lineLengths)) / 2, IMPEDANCE_START_Y, (centreX - lineLengths) + ((centreX - impedanceWidth / 2)  - (centreX - lineLengths)) / 2 - currentArraySize, IMPEDANCE_START_Y - currentArraySize);

  // text labels
  fill(200);
  textAlign(CENTER);
  text("Z = R + jX", centreX, IMPEDANCE_START_Y - 25);
  fill(VaColor);
  text("Va", centreX - lineLengths, IMPEDANCE_START_Y - 25);
  fill(VbColor);
  text("Vb", centreX + lineLengths, IMPEDANCE_START_Y - 25);
  textAlign(LEFT);
  fill(IColor);
  text("i", (centreX - lineLengths) + ((centreX - impedanceWidth / 2)  - (centreX - lineLengths)) / 2 - currentArraySize, IMPEDANCE_START_Y + currentArraySize + 15);
}

void drawImpedancePlot() {
  // draw axes
  stroke(50);
  strokeWeight(2);
  line(PHASORS_START_X - ONE_PU_HEIGHT / 4, IMPEDANCE_START_Y, PHASORS_START_X + 2 * ONE_PU_HEIGHT, IMPEDANCE_START_Y);
  line(PHASORS_START_X, IMPEDANCE_START_Y - 1.5 * ONE_PU_HEIGHT, PHASORS_START_X, IMPEDANCE_START_Y + 1.5 * ONE_PU_HEIGHT);

  text("R", PHASORS_START_X + 2 * ONE_PU_HEIGHT + 1, IMPEDANCE_START_Y + 11);
  text("+X", PHASORS_START_X, IMPEDANCE_START_Y - 1.5 * ONE_PU_HEIGHT + 11);
  text("-X", PHASORS_START_X, IMPEDANCE_START_Y + 1.5 * ONE_PU_HEIGHT + 11);

  ellipseMode(CENTER);
  ellipse(PHASORS_START_X + R * ONE_PU_HEIGHT, IMPEDANCE_START_Y - X * ONE_PU_HEIGHT, 7, 7);
}

void drawPhasors(float m, float ang, int x, int y, color c) {
  // draw axes
  stroke(50);
  strokeWeight(2);
  line(x - ONE_PU_HEIGHT, y, x + ONE_PU_HEIGHT, y);
  line(x, y - ONE_PU_HEIGHT, x, y + ONE_PU_HEIGHT);

  float mDraw = m * ONE_PU_HEIGHT;
  ang = radians(ang);

  for (int j = 0; j < 3; j++) {
    // shift phases B and C by 120deg (-/+ respectively)
    if (j == 1) {
      ang -= radians(120);
    }
    else if (j == 2) {
      ang += radians(240);
    }

    strokeWeight(OFFSET_ORIGINAL);
    stroke(c, 255 - j * 100);
    line(x, y, x + mDraw*cos(ang), y - mDraw*sin(ang));

    if (mouseIsOverLine(x, y, x + mDraw*cos(ang), y - mDraw*sin(ang))) {
      fill(210);
      text(nf(m, 1, 2) + "∠ " + nf(degrees(ang), 1) + "°", x + mDraw*cos(ang), y - mDraw*sin(ang));
    }
  }
}

float mean(float[] values) {
  float total = 0.0;

  for (int index = 0; index < ITERATIONS; index++) {
    total = total + values[index];
  }

  return total / ITERATIONS;
}

float rms(float[] values) {
  float squareTotal = 0.0;
  int PERIOD_ITERATIONS = 1;

  if (totalTime > (1 / f)) {
    PERIOD_ITERATIONS = int(1 / (f * Ts));
  }
  else {
    PERIOD_ITERATIONS = ITERATIONS;
  }

  for (int index = 0; index < PERIOD_ITERATIONS; index++) {
    squareTotal = squareTotal + (values[index] * values[index]);
  }

  return sqrt(squareTotal / float(PERIOD_ITERATIONS));
}

void plot(float[] data, int startY, int c, String name, boolean showMean, boolean showRMS) {
  // show name
  fill(c);
  //textFont(fontBold, 16);
  text(name, PHASORS_START_X - 90 + (name.equals("Q") ? 40 : 0), startY);
  //textFont(font, 16);

  // draw x-axis
  stroke(50);
  strokeWeight(2);
  line(WAVEFORMS_START_X, startY, WAVEFORMS_END_X, startY);

  // draw waveform
  strokeWeight(2);
  stroke(c, 255);
  noFill();
  beginShape(POLYGON);

  for (int t = 0; t < ITERATIONS; t++) {
    vertex(WAVEFORMS_START_X + (float(t) * float(WAVEFORMS_END_X - WAVEFORMS_START_X)) / ITERATIONS, startY - (data[t] * ONE_PU_HEIGHT));
  }

  endShape();

  // add labels
  if (showRMS == true) {
    text("RMS = " + nf(rms(data), 1, 3), WAVEFORMS_END_X + 10, startY);
  }

  if (showMean == true) {
    text("mean = " + nf(mean(data), 1, 3), WAVEFORMS_END_X + 10, startY/* + 15*/);
  }
}

void setup() {
  size(MAX_WIDTH, MAX_HEIGHT/*, OPENGL*/);
  background(0);
  randomSeed(0);
  //strokeWeight(2);
  strokeCap(ROUND);

  font = loadFont("Courier New Bold", 16);
  //font = loadFont("monospace", 16);
  textFont(font, 16);
  textLeading(10);
  textAlign(LEFT);
  smooth();
  colorMode(HSB);
  frameRate(30);

  reset();	// implemented in JavaScript
}

// returns phase (in degrees) of phase of sinusoidal data x[]
float getPhase(float[] x) {
    float wavepoint = x[0] / (max(rms(x), 0.001) * sqrt(2));	// max() prevents divide by zero

    // ensure valid range before calling asin()
    if (wavepoint < -1.0) {
      wavepoint = -1.0;
    }
    else if (wavepoint > 1.0) {
      wavepoint = 1.0;
    }

    phase = degrees(asin(wavepoint));

    //
    float dx = (x[1] - x[0]);
    if (dx < 0) {
      phase = 180 - phase;
    }
    
    return phase;
}

void draw() {
  VaMag = Processing.data.VaMag;
  VbMag = Processing.data.VbMag;
  VaPhase = Processing.data.VaPhase;
  VbPhase = Processing.data.VbPhase;
  R = Processing.data.R;
  X = Processing.data.X;
  totalTime = Processing.data.t;
  f = Processing.data.f;

  if (Processing.data.change == true) {
    background(0);

    Processing.data.change = false;

    updateLabels();   // implemented in JavaScript

    ITERATIONS = int(totalTime / Ts);
    w = 2 * PI * f;
    L = X / w;
    C = -1 / (X * w);
    Zmag = sqrt(R*R + X*X);

    // compute waveforms at each time-step
    for (int t = 0; t < ITERATIONS; t++) {
      Va[t] = VaMag*sin(w*(float(t) * Ts) + (VaPhase * PI / 180.0));
      Vb[t] = VbMag*sin(w*(float(t) * Ts) + (VbPhase * PI / 180.0));

      dV[t] = Va[t] - Vb[t];

      if (X > 0) {
        i[t] = (1 / (L*L * w*w + R*R)) * (-L * VaMag * w * cos(w*(float(t) * Ts) + (VaPhase * PI / 180.0)) + R*VaMag*sin(w*(float(t) * Ts) + (VaPhase * PI / 180.0)) + L * VbMag * w * cos(w*(float(t) * Ts) + (VbPhase * PI / 180.0)) - R*VbMag*sin(w*(float(t) * Ts) + (VbPhase * PI / 180.0)));
      }
      else if (X < 0) {
        i[t] = (C*R*VaMag*w*cos((VaPhase * PI / 180.0) + (float(t) * Ts)*w) - C*R*VbMag*w*cos((VbPhase * PI / 180.0) + (float(t) * Ts)*w) + C*C*R*R*VaMag*w*w*sin((VaPhase * PI / 180.0) + (float(t) * Ts)*w) - C*C*R*R*VbMag*w*w*sin((VbPhase * PI / 180.0) + (float(t) * Ts)*w)) / (R*(C*C*R*R*w*w + 1));
      }
      else {
        i[t] = dV[t] / R;
      }

      P[t] = dV[t] * i[t];			// or P[t] = pow(i[t], 2) * R;
      //P[t] = i[t] * i[t] * R;
      //P[t] = dV[t] * dV[t] / Zmag;
      Q[t] = i[t] * i[t] * X;      //Q[t] = pow(dV[t], 2) / X;
      //Q[t] = dV[t] * dV[t] * X;      //Q[t] = pow(dV[t], 2) / X;
    }

    meanP = mean(P);
    meanQ = mean(Q);
    rmsI = rms(i);
    rmsdV = rms(dV);

    float S = sqrt(meanP * meanP + meanQ * meanQ);
    float pf = meanP / S;
    float IPhase = getPhase(i);
    float dVPhase = getPhase(dV);
    
    if (Processing.data.PQCombined == false) {
      for (int t = 0; t < ITERATIONS; t++) {
        P[t] = rmsdV * rmsI * pf * (1 + cos(-PI + 2 * (dVPhase * PI / 180.0) + 2 * w * (float(t) * Ts)));
        Q[t] = rmsdV * rmsI * cos(asin(pf)) * sin(-PI + 2 * (dVPhase * PI / 180.0) + 2 * w * (float(t) * Ts));
        if (X < 0) {
          Q[t] = -1 * Q[t];
        }
      }
    }

    drawImpedancePlot();
    drawPowerSystem();

    plot(Va, VA_START, VaColor, "Va", false, true);
    plot(Vb, VB_START, VbColor, "Vb", false, true);
    plot(i, I_START, IColor, "i", false, true);
    //plot(P, P_START, PColor, "P + Q", true);
    
    if (Processing.data.PQCombined === false) {
      plot(P, P_START, PColor, "P", false, false);
      plot(Q, P_START, QColor, "Q", false, false);
    }
    else {
      plot(P, P_START, PColor, "P + Q", true, false);
    }

    colorMode(RGB);
    fill(255, 255, 255, 200);
    colorMode(HSB);
    text("P = " + nf(rmsdV * rmsI * pf, 1, 3) + " W", WAVEFORMS_END_X + 10, P_START + 20);
    text("Q = " + nf(rmsdV * rmsI * cos(asin(pf)), 1, 3) + " var", WAVEFORMS_END_X + 10, P_START + 35);
    text("|S| = " + nf(S, 1, 3) + " VA", WAVEFORMS_END_X + 10, P_START + 50);
    text("p.f. = " + nf(pf, 1, 3) + " " + (X < 0 ? "lead" : (X > 0 ? "lag" : "")), WAVEFORMS_END_X + 10, P_START + 65);

    drawPhasors(VaMag, VaPhase, PHASORS_START_X, VA_START, VaColor);
    drawPhasors(VbMag, VbPhase, PHASORS_START_X, VB_START, VbColor);
    drawPhasors(rmsI * sqrt(2), IPhase, PHASORS_START_X, I_START, IColor);

    if (X > 0) {
      text("L = " + nf(L * 1000, 1, 3) + "mH", WAVEFORMS_END_X + 10, P_START + 80);
    }
    else if (X < 0) {
      text("C = " + nf(C * 1000, 1, 3) + "mF", WAVEFORMS_END_X + 10, P_START + 80);
    }
  }
}

