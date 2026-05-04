(function () {
  'use strict';

  const dom = {};
  const state = {
    prepared: false,
    running: false,
    aborted: false,
    breakActive: false,
    stageWidth: 0,
    stageHeight: 0,
    renderWidth: 0,
    renderHeight: 0,
    displayScaleX: 1,
    displayScaleY: 1,
    sessionSeed: 0,
    sessionStamp: null,
    settings: null,
    trialSpecs: [],
    results: [],
    activeTrial: null,
    activeKeyResolver: null,
    activeResponseTimeoutId: null,
    autoDownloadDone: false,
    autoPrepareTimerId: null,
    serverShutdownTimerId: null,
    instructionDemoSource: null,
    instructionDemoSpec: null,
    quest: null,
    questHistory: [],
  };

  const DEFAULTS = {
    trialCount: 40,
    fixMinMs: 800,
    fixMaxMs: 1400,
    responseTimeoutMs: 5000,
    bgContrast: 0.22,
    targetContrast: 0.28,
    rtCriterionMs: 2700,
    noiseDiskHeightRatio: 0.9,
    responseKey: 'Space',
    flashMs: 700,
    itiMs: 450,
    bgWidthDeg: 15,
    gaborSF: 4,
    gaborCyc: 2,
    gaborWidthDeg: 0.3,
    gaborOrientDeg: -45,
    eccMinDeg: 2,
    eccMaxDeg: 7,
    learnP: 0.5,
    testP: 0.4,
    connectTNum: 50,
    restInterval: 20,
    bgFlashDurationMs: 400,   // 背景闪烁总时长（2次）
  };

  document.addEventListener('DOMContentLoaded', init);

  function init() {
    cacheDom();
    wireEvents();
    document.getElementById('appShell').style.display = 'none';
    document.getElementById('startOverlay').style.display = 'flex';
    // 预绘制指导图
    drawFixPointDemo();
    drawFixMiniDemo();
    drawTargetDemo();
    resizeStageCanvas();
    drawIdleScreen();
    setStatus('等待准备。');
    setSummary('暂无结果。');
    initSubjectCounter();
    updateButtons();
      armAutoFullscreen();
      requestFullscreenSafe();
  }

  function cacheDom() {
    dom.overlay = document.getElementById('startOverlay');
    dom.infoSection = document.getElementById('infoSection');
    dom.instrSection = document.getElementById('instrSection');
    dom.overlayNextBtn = document.getElementById('overlayNextBtn');
    dom.instrBackBtn = document.getElementById('instrBackBtn');
    dom.instrStartBtn = document.getElementById('instrStartBtn');
    dom.finishSection = document.getElementById('finishSection');
    dom.finishCloseBtn = document.getElementById('finishCloseBtn');
    dom.finishDownloadJsonBtn = document.getElementById('finishDownloadJsonBtn');
    dom.finishDownloadCsvBtn = document.getElementById('finishDownloadCsvBtn');
    dom.finishDownloadTxtBtn = document.getElementById('finishDownloadTxtBtn');
    dom.finishCopyJsonBtn = document.getElementById('finishCopyJsonBtn');

    dom.subjectId = document.getElementById('overlaySubjectId');
    dom.subjectName = document.getElementById('overlaySubjectName');
    dom.subjectGender = document.getElementById('overlaySubjectGender');
    dom.organization = document.getElementById('overlayOrganization');
    dom.sessionId = document.getElementById('overlaySessionId');
    dom.restInterval = document.getElementById('overlayRestInterval');
    dom.trialCount = document.getElementById('overlayTrialCount');
    dom.fixMinMs = document.getElementById('overlayFixMinMs');
    dom.fixMaxMs = document.getElementById('overlayFixMaxMs');
    dom.responseTimeoutMs = document.getElementById('overlayResponseTimeoutMs');
    dom.rtCriterionMs = document.getElementById('overlayRtCriterionMs');
    dom.bgContrast = document.getElementById('overlayBgContrast');
    dom.targetContrast = document.getElementById('overlayTargetContrast');
    dom.noiseDiskHeightRatio = document.getElementById('overlayNoiseDiskHeightRatio');
    dom.responseKey = document.getElementById('overlayResponseKey');

    dom.sidebar = document.getElementById('sidebar');
    dom.stagePanel = document.getElementById('stagePanel');
    dom.stageFrame = document.getElementById('stageFrame');
    dom.canvas = document.getElementById('experimentCanvas');
    dom.ctx = dom.canvas.getContext('2d');
    dom.statusText = document.getElementById('statusText');
    dom.summaryText = document.getElementById('summaryText');
    dom.prepProgressFill = document.getElementById('prepProgressFill');
    dom.startBtn = document.getElementById('startBtn');

    dom.instrCanvas1 = document.getElementById('instrCanvas1');
    dom.instrCanvas1Solid = document.getElementById('instrCanvas1Solid');
    dom.instrCanvas2 = document.getElementById('instrCanvas2');
    dom.instrFixMiniHollow = document.getElementById('instrFixMiniHollow');
    dom.instrFixMiniSolid = document.getElementById('instrFixMiniSolid');
    dom.instrTargetMini = document.getElementById('instrTargetMini');
  }

  function wireEvents() {
    dom.overlayNextBtn.addEventListener('click', () => {
      dom.infoSection.style.display = 'none';
      dom.instrSection.style.display = 'block';
    });

    dom.instrBackBtn.addEventListener('click', () => {
      dom.instrSection.style.display = 'none';
      dom.infoSection.style.display = 'block';
    });

    dom.instrStartBtn.addEventListener('click', () => {
      dom.overlay.style.display = 'none';
      document.getElementById('appShell').style.display = 'block';
      requestFullscreenSafe();
      scheduleAutoPrepare(0);
      setTimeout(() => {
        if (state.prepared) dom.startBtn.click();
        else {
          const check = () => {
            if (state.prepared) dom.startBtn.click();
            else if (!state.aborted) requestAnimationFrame(check);
          };
          requestAnimationFrame(check);
        }
      }, 150);
    });

    if (dom.finishCloseBtn) {
      dom.finishCloseBtn.addEventListener('click', () => {
        dom.overlay.style.display = 'none';
      });
    }
    if (dom.finishDownloadJsonBtn) {
      dom.finishDownloadJsonBtn.addEventListener('click', () => {
        downloadResults('json');
        setStatus('已触发 JSON 下载。');
      });
    }
    if (dom.finishDownloadCsvBtn) {
      dom.finishDownloadCsvBtn.addEventListener('click', () => {
        downloadResults('csv');
        setStatus('已触发 CSV 下载。');
      });
    }
    if (dom.finishDownloadTxtBtn) {
      dom.finishDownloadTxtBtn.addEventListener('click', () => {
        downloadResults('txt');
        setStatus('已触发 TXT 下载。');
      });
    }
    if (dom.finishCopyJsonBtn) {
      dom.finishCopyJsonBtn.addEventListener('click', () => {
        copyJsonToClipboard();
      });
    }

    dom.startBtn.addEventListener('click', () => {
      startExperiment().catch(handleFatalError);
    });

    const paramInputs = [
      dom.trialCount, dom.fixMinMs, dom.fixMaxMs,
      dom.responseTimeoutMs, dom.rtCriterionMs,
      dom.bgContrast, dom.targetContrast, dom.noiseDiskHeightRatio,
      dom.responseKey, dom.sessionId, dom.restInterval
    ];
    paramInputs.forEach(el => {
      if (el) {
        el.addEventListener('input', () => scheduleAutoPrepare(180));
        el.addEventListener('change', () => scheduleAutoPrepare(180));
      }
    });

    document.addEventListener('keydown', handleKeyDown, true);
    window.addEventListener('resize', () => {
      resizeStageCanvas();
      if (!state.running && !state.breakActive) drawIdleScreen();
    });
  }

  function readSettings() {
    const trialCount = clampInt(parseInt(dom.trialCount.value), 1, 500, DEFAULTS.trialCount);
    const fixMinMs = clampInt(parseInt(dom.fixMinMs.value), 100, 20000, DEFAULTS.fixMinMs);
    const fixMaxMs = clampInt(parseInt(dom.fixMaxMs.value), 100, 20000, DEFAULTS.fixMaxMs);
    const responseTimeoutMs = clampInt(parseInt(dom.responseTimeoutMs.value), 500, 60000, DEFAULTS.responseTimeoutMs);
    const rtCriterionMs = clampInt(parseInt(dom.rtCriterionMs.value), 100, 20000, DEFAULTS.rtCriterionMs);
    const bgContrast = clampFloat(parseFloat(dom.bgContrast.value), 0, 1, DEFAULTS.bgContrast);
    const targetContrast = clampFloat(parseFloat(dom.targetContrast.value), 0, 1, DEFAULTS.targetContrast);
    const noiseDiskHeightRatio = clampFloat(parseFloat(dom.noiseDiskHeightRatio.value), 0.6, 1.0, DEFAULTS.noiseDiskHeightRatio);
    const restInterval = clampInt(parseInt(dom.restInterval.value), 1, 500, DEFAULTS.restInterval);

    return {
      subjectId: String(dom.subjectId.value || '1').trim(),
      subjectName: String(dom.subjectName.value || '').trim(),
      subjectGender: String(dom.subjectGender.value || '').trim(),
      organization: String(dom.organization.value || '').trim(),
      sessionId: String(dom.sessionId.value || '1').trim(),
      trialCount, fixMinMs: Math.min(fixMinMs, fixMaxMs), fixMaxMs: Math.max(fixMinMs, fixMaxMs),
      responseTimeoutMs, rtCriterionMs, bgContrast, targetContrast, noiseDiskHeightRatio,
      responseKey: String(dom.responseKey.value || DEFAULTS.responseKey),
      flashMs: DEFAULTS.flashMs, itiMs: DEFAULTS.itiMs,
      bgWidthDeg: DEFAULTS.bgWidthDeg, gaborSF: DEFAULTS.gaborSF,
        gaborCyc: DEFAULTS.gaborCyc, gaborWidthDeg: DEFAULTS.gaborWidthDeg, gaborOrientDeg: DEFAULTS.gaborOrientDeg,
      eccMinDeg: DEFAULTS.eccMinDeg, eccMaxDeg: DEFAULTS.eccMaxDeg,
      learnP: DEFAULTS.learnP, testP: DEFAULTS.testP, connectTNum: DEFAULTS.connectTNum,
      restInterval, bgFlashDurationMs: DEFAULTS.bgFlashDurationMs,
      timeZone: Intl.DateTimeFormat().resolvedOptions().timeZone || '',
      locale: navigator.language || '',
      userAgent: navigator.userAgent || '',
      screenWidth: window.screen.width, screenHeight: window.screen.height,
      devicePixelRatio: window.devicePixelRatio || 1,
      pageWidth: window.innerWidth, pageHeight: window.innerHeight,
    };
  }

  async function prepareExperiment() {
    state.settings = readSettings();
    state.sessionSeed = createSessionSeed(state.settings);
    state.results = [];
    state.trialSpecs = [];
    state.questHistory = [];
    state.quest = null;
    state.prepared = false;
    state.running = false;
    state.aborted = false;
    state.autoDownloadDone = false;
    state.breakActive = false;
    state.instructionDemoSource = null;
    state.instructionDemoSpec = null;
    cancelServerShutdown();

    resizeStageCanvas();
    const available = getStageCssSize();
    const dpr = window.devicePixelRatio || 1;
    state.renderWidth = Math.max(1, Math.round(available.width * dpr));
    state.renderHeight = Math.max(1, Math.round(available.height * dpr));
    state.displayScaleX = available.width / state.renderWidth;
    state.displayScaleY = available.height / state.renderHeight;

    appendStatus(`准备刺激中... Canvas: ${available.width}×${available.height}`);
    setProgress(0);
    state.trialSpecs = buildTrialSpecs(state.settings, state.sessionSeed, state.renderWidth, state.renderHeight);
    state.quest = initializeQuestState(state.settings, state.trialSpecs.length);
    setProgress(1);
    state.prepared = !state.aborted;
    if (state.prepared) {
      setStatus('刺激准备完成，可以开始实验。');
      dom.startBtn.disabled = false;
      drawIdleScreen();
    }
    updateButtons();
  }

  async function startExperiment() {
    await prepareExperiment();
    if (!state.prepared) return;

    document.body.classList.add('running');
    document.body.style.cursor = 'none';  // 隐藏鼠标
    if (dom.overlay) dom.overlay.style.display = 'none';
    await nextAnimationFrame();
    resizeStageCanvas();
    rebuildGeometryForCurrentStage();

    state.running = true;
    state.aborted = false;
    state.breakActive = false;
    updateButtons();
    state.sessionStamp = makeStamp(performance.now());
    appendStatus(`实验开始。\nSession start: ${state.sessionStamp.localCompact}`);

    for (let i = 0; i < state.trialSpecs.length; i++) {
      if (state.aborted) break;
      const spec = state.trialSpecs[i];
      const questPlan = computeQuestTrialPlan(state.quest, state.questHistory, i + 1, state.trialSpecs.length);
      spec.targetContrast = questPlan.targetContrast;
      spec.quest = questPlan;
      const result = await runSingleTrial(i, spec);
      const judge = getQuestJudgeFromTrial(result, state.settings.rtCriterionMs);
      result.questJudge = judge;
      result.questRtCriterionMs = state.settings.rtCriterionMs;
      updateQuestPosterior(state.quest.q, questPlan.tTest, judge);
      state.questHistory.push(judge);
      state.results.push(result);
      updateLiveSummary(result);
      setProgress((i + 1) / state.trialSpecs.length);
      setStatus(`运行中... (${i + 1}/${state.trialSpecs.length})`);

      const restInterval = state.settings.restInterval;
      if ((i + 1) % restInterval === 0 && (i + 1) < state.trialSpecs.length) {
        await takeBreak();
        if (state.aborted) break;
      }
      await waitMs(state.settings.itiMs);
    }

    state.running = false;
    document.body.classList.remove('running');
    document.body.style.cursor = '';
    await exitFullscreenSafe();
    drawIdleScreen();
    if (!state.aborted) {
      setStatus('实验完成，数据已自动下载。');
      setSummary(buildSummaryText());
      autoDownloadAllResults();
      incrementSubjectCounter();
      state.prepared = false;
      showFinishOverlay();
      scheduleServerShutdown(4500);
    } else {
      setStatus('实验已中止。');
    }
    updateButtons();
  }

  // ---------- 重构的 trial 流程 ----------
  async function runSingleTrial(trialIndex, spec) {
    const trialCanvas = buildTrialCanvas(spec, state.renderWidth, state.renderHeight);
    const result = {
      trialIndex: trialIndex + 1,
      seed: spec.seed, hueDeg: round3(spec.hueDeg), saturation: round3(spec.saturation),
      bgContrast: round3(spec.bgContrast), targetContrast: round3(spec.targetContrast),
      fixationHoldMs: Math.round(spec.fixationHoldMs),
      responseTimeoutMs: spec.responseTimeoutMs,
      renderWidth: state.renderWidth, renderHeight: state.renderHeight,
      stageWidth: state.stageWidth, stageHeight: state.stageHeight,
      targetXpx: round3(spec.targetX), targetYpx: round3(spec.targetY),
      targetXnorm: round4(spec.targetX / state.renderWidth),
      targetYnorm: round4(spec.targetY / state.renderHeight),
      targetRadiusPx: round3(spec.targetRadiusPx), targetSigmaPx: round3(spec.targetSigmaPx),
      targetOrientationDeg: round3(spec.targetOriDeg),
      carrierCyclesPerPx: round6(spec.carrierCyclesPerPx),
      bgRadiusPx: round3(spec.bgRadiusPx),
      fixationPress: null,      // 新增
      fixationSolid: null,      // 新增
      stimulus: null,
      response: null,
      backgroundFlash: null,    // 新增
      flash: null,
      trialEnd: null,
      respondAfterFixation: false,
      responded: false,
      timeout: false,
    };

    state.activeTrial = result;

    // 1. 绘制空心注视点，等待按键
    drawFixationPhase('hollow');
    const fixationPress = await waitForResponse(null); // 无超时，必须按键后继续
    result.fixationPress = makeStamp(fixationPress.perfNow);
    result.respondAfterFixation = true;

    // 2. 注视点变实心，并开始随机等待
    drawFixationPhase('solid');
    const fixationSolidStamp = performance.now();
    result.fixationSolid = makeStamp(fixationSolidStamp);
    await waitMs(spec.fixationHoldMs);

    // 3. 显示刺激
    drawTrialFrame(trialCanvas, spec, 'stimulus');
    const stimStamp = await nextAnimationFrame();
    result.stimulus = makeStamp(stimStamp);
    appendTrialLog(result, `Stimulus on ${result.stimulus.localCompact}`);

    // 4. 等待按键响应（或超时）
    const response = await waitForResponse(spec.responseTimeoutMs);
    result.response = response.timestamp;
    result.responded = response.responded;
    result.timeout = response.timedOut;
    result.responseKey = response.keyCode || '';
    result.responseKeyLabel = response.keyLabel || '';
    result.responseRTMs = response.responded ? round3(response.perfNow - stimStamp) : null;
      result.questTargetP = spec.quest ? round4(spec.quest.pThreshold) : null;
      result.questGrain = spec.quest ? round4(spec.quest.grain) : null;
      result.questBoostLog10 = spec.quest ? round4(spec.quest.boostLog10) : null;
      result.questIntensityLog10 = spec.quest ? round4(spec.quest.tTest) : null;

    // 5. 如果按键，先执行背景闪烁反馈，再执行光圈闪烁
    if (response.responded) {
      const bgFlashEnd = await performBackgroundFlash(trialCanvas, spec);
      result.backgroundFlash = makeStamp(bgFlashEnd);
      const flashEnd = await flashTarget(trialCanvas, spec);
      result.flash = makeStamp(flashEnd);
      result.flashed = true;
    }

    result.trialEnd = makeStamp(performance.now());
    state.activeTrial = null;
    return result;
  }

  // 绘制注视点（空心或实心）
  function drawFixationPhase(type) {
    const ctx = dom.ctx;
    const cssW = state.stageWidth;
    const cssH = state.stageHeight;
    ctx.save();
    ctx.setTransform(window.devicePixelRatio || 1, 0, 0, window.devicePixelRatio || 1, 0, 0);
    ctx.clearRect(0, 0, cssW, cssH);
    ctx.fillStyle = '#808080';
    ctx.fillRect(0, 0, cssW, cssH);

    const cx = cssW / 2, cy = cssH / 2;
    const radius = Math.max(6, Math.round(Math.min(cssW, cssH) * 0.008));
    ctx.beginPath();
    ctx.arc(cx, cy, radius, 0, Math.PI * 2);
    if (type === 'hollow') {
      ctx.strokeStyle = '#000';
      ctx.lineWidth = 3;
      ctx.stroke();
    } else {
      ctx.fillStyle = '#000';
      ctx.fill();
    }
    ctx.restore();
  }

  // 背景闪烁反馈（整个背景变浅/恢复 2 次）
  function performBackgroundFlash(trialCanvas, spec) {
    return new Promise(resolve => {
      const ctx = dom.ctx;
      const cssW = state.stageWidth;
      const cssH = state.stageHeight;
      const flashDuration = state.settings.bgFlashDurationMs;
      const halfPeriod = flashDuration / 4; // 两次闪烁 = 4个阶段
      const startTime = performance.now();
      const cx = cssW / 2;
      const cy = cssH / 2;
      const scale = Math.min(state.displayScaleX, state.displayScaleY);
      const flashRadius = Math.max(8, spec.bgRadiusPx * scale);

      function drawFlash(now) {
        const elapsed = now - startTime;
        const period = Math.floor(elapsed / halfPeriod);
        const showOverlay = period % 2 === 0; // 0,2 -> 遮罩；1,3 -> 正常
        ctx.save();
        ctx.setTransform(window.devicePixelRatio || 1, 0, 0, window.devicePixelRatio || 1, 0, 0);
        ctx.clearRect(0, 0, cssW, cssH);
        ctx.drawImage(trialCanvas, 0, 0, cssW, cssH);
        if (showOverlay) {
          ctx.save();
          ctx.beginPath();
          ctx.arc(cx, cy, flashRadius, 0, Math.PI * 2);
          ctx.clip();
          ctx.fillStyle = 'rgba(255, 255, 255, 0.3)';
          ctx.fillRect(cx - flashRadius - 2, cy - flashRadius - 2, flashRadius * 2 + 4, flashRadius * 2 + 4);
          ctx.restore();
        }
        ctx.restore();

        if (elapsed < flashDuration) {
          requestAnimationFrame(drawFlash);
        } else {
          // 最后恢复无遮罩画面
          ctx.save();
          ctx.setTransform(window.devicePixelRatio || 1, 0, 0, window.devicePixelRatio || 1, 0, 0);
          ctx.clearRect(0, 0, cssW, cssH);
          ctx.drawImage(trialCanvas, 0, 0, cssW, cssH);
          ctx.restore();
          resolve(now);
        }
      }
      requestAnimationFrame(drawFlash);
    });
  }

  async function takeBreak() {
    state.breakActive = true;
    setStatus('休息时间，按回车键继续...');
    drawBreakScreen();
    return new Promise((resolve) => {
      const onKey = (e) => {
        if (e.code === 'Enter') {
          e.preventDefault();
          document.removeEventListener('keydown', onKey, true);
          state.breakActive = false;
          resolve();
        } else if (e.code === 'Escape') {
          e.preventDefault();
          document.removeEventListener('keydown', onKey, true);
          abortExperiment('用户按下 Escape 中止。');
          resolve();
        }
      };
      document.addEventListener('keydown', onKey, true);
    });
  }

  function drawBreakScreen() {
    resizeStageCanvas();
    const ctx = dom.ctx;
    const cssW = state.stageWidth;
    const cssH = state.stageHeight;
    ctx.save();
    ctx.setTransform(window.devicePixelRatio || 1, 0, 0, window.devicePixelRatio || 1, 0, 0);
    ctx.clearRect(0, 0, cssW, cssH);
    ctx.fillStyle = 'rgb(128,128,128)';
    ctx.fillRect(0, 0, cssW, cssH);
    ctx.fillStyle = 'rgba(0,0,0,0.8)';
    ctx.font = 'bold 28px Georgia, serif';
    ctx.textAlign = 'center';
    ctx.fillText('休息一下', cssW/2, cssH/2 - 30);
    ctx.font = '20px "Trebuchet MS", sans-serif';
    ctx.fillText('按回车键继续', cssW/2, cssH/2 + 20);
    ctx.restore();
  }

  // 指导语示例图
  function drawFixPointDemo() {
    const drawFixDemo = (canvas, solid, label) => {
      if (!canvas) return;
      const ctx = canvas.getContext('2d');
      const w = canvas.width;
      const h = canvas.height;
      ctx.clearRect(0, 0, w, h);
      ctx.fillStyle = '#808080';
      ctx.fillRect(0, 0, w, h);
      ctx.beginPath();
      ctx.arc(w / 2, h / 2, 12, 0, Math.PI * 2);
      if (solid) {
        ctx.fillStyle = '#000';
        ctx.fill();
      } else {
        ctx.strokeStyle = '#000';
        ctx.lineWidth = 3;
        ctx.stroke();
      }
      ctx.font = '10px sans-serif';
      ctx.fillStyle = '#000';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'alphabetic';
      ctx.fillText(label, w / 2, h - 10);
    };

    drawFixDemo(dom.instrCanvas1Solid, true, '等待开始');
    drawFixDemo(dom.instrCanvas1, false, '按空格键');
  }

  function drawFixMiniDemo() {
    const drawMini = (canvas, solid) => {
      if (!canvas) return;
      const ctx = canvas.getContext('2d');
      const w = canvas.width;
      const h = canvas.height;
      ctx.fillStyle = '#808080';
      ctx.fillRect(0, 0, w, h);
      ctx.beginPath();
      ctx.arc(w / 2, h / 2, 8, 0, Math.PI * 2);
      if (solid) {
        ctx.fillStyle = '#000';
        ctx.fill();
      } else {
        ctx.strokeStyle = '#000';
        ctx.lineWidth = 2;
        ctx.stroke();
      }
    };
    drawMini(dom.instrFixMiniHollow, false);
    drawMini(dom.instrFixMiniSolid, true);
  }

  function drawTargetDemo() {
    const canvas = dom.instrCanvas2;
    if (!canvas) return;
    const source = buildInstructionDemoSource();
    state.instructionDemoSource = source;
    drawScaledCanvas(canvas, source, 0, 0, source.width, source.height);
    drawTargetMiniDemo();
  }

  function drawTargetMiniDemo() {
    const canvas = dom.instrTargetMini;
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    const w = canvas.width;
    const h = canvas.height;
    const off = document.createElement('canvas');
    off.width = w;
    off.height = h;
    const offCtx = off.getContext('2d', { willReadFrequently: true });

    // 独立绘制 mini target，避免从大背景裁剪导致目标过小
    const cx = w / 2;
    const cy = h / 2;
    const sigma = Math.min(w, h) * 0.27;
    const contrast = 0.6;
    const freq = 0.065; 
    const ori = -45;
    const theta = degToRad(ori);
    const cosT = Math.cos(theta);
    const sinT = Math.sin(theta);
    const sigmaSq2 = 2 * sigma * sigma;
    const img = offCtx.createImageData(w, h);

    for (let y = 0; y < h; y += 1) {
      for (let x = 0; x < w; x += 1) {
        const tx = x - cx;
        const ty = y - cy;
        const xr = tx * cosT + ty * sinT;
        const env = Math.exp(-((tx * tx) + (ty * ty)) / sigmaSq2);
        const carrier = Math.cos(2 * Math.PI * freq * xr);
        const val = clamp01(0.5 + contrast * env * carrier);
        const idx = (y * w + x) * 4;
        const g = clamp255(val * 255);
        img.data[idx] = g;
        img.data[idx + 1] = g;
        img.data[idx + 2] = g;
        img.data[idx + 3] = 255;
      }
    }
    offCtx.putImageData(img, 0, 0);

    ctx.clearRect(0, 0, w, h);
    ctx.save();
    const r = Math.min(w, h) * 0.48;
    ctx.beginPath();
    ctx.arc(cx, cy, r, 0, Math.PI * 2);
    ctx.clip();
    ctx.drawImage(off, 0, 0);
    ctx.restore();
  }


  function buildInstructionDemoSource() {
    const previewCanvas = dom.instrCanvas2;
    const baseWidth = Math.max(960, Math.round((previewCanvas ? previewCanvas.width : 360) * 3));
    const baseHeight = Math.max(640, Math.round(baseWidth * 2 / 3));
    const demoSettings = {
      ...DEFAULTS,
      trialCount: 1,
      bgContrast: 0.22,
      targetContrast: 0.0,
      noiseDiskHeightRatio: 0.9,
    };
    const spec = buildTrialSpecs(demoSettings, 999, baseWidth, baseHeight)[0];
    spec.hueDeg = 80;        // 色调（范围 0~360）
    state.instructionDemoSpec = spec;
    return buildTrialCanvas(spec, baseWidth, baseHeight);
  }

  function buildInstructionDemoCrop(source, destWidth, destHeight) {
    const spec = state.instructionDemoSpec;
    const targetX = spec && Number.isFinite(spec.targetX) ? spec.targetX : source.width * 0.68;
    const targetY = spec && Number.isFinite(spec.targetY) ? spec.targetY : source.height * 0.42;
    const cropWidth = Math.min(source.width, Math.max(Math.round(destWidth * 1.25), 180));
    const cropHeight = Math.min(source.height, Math.max(Math.round(destHeight * 1.25), 120));
    const x = clampFloat(targetX - cropWidth / 2, 0, Math.max(0, source.width - cropWidth), targetX - cropWidth / 2);
    const y = clampFloat(targetY - cropHeight / 2, 0, Math.max(0, source.height - cropHeight), targetY - cropHeight / 2);
    return { x, y, width: cropWidth, height: cropHeight };
  }

  function drawScaledCanvas(targetCanvas, sourceCanvas, sx, sy, sw, sh) {
    const ctx = targetCanvas.getContext('2d');
    const w = targetCanvas.width;
    const h = targetCanvas.height;
    ctx.save();
    ctx.clearRect(0, 0, w, h);
    ctx.imageSmoothingEnabled = true;
    ctx.imageSmoothingQuality = 'high';
    ctx.drawImage(sourceCanvas, sx, sy, sw, sh, 0, 0, w, h);
    ctx.restore();
  }


  function drawTrialFrame(trialCanvas, spec, phase) {
    const ctx = dom.ctx;
    const cssW = state.stageWidth;
    const cssH = state.stageHeight;
    ctx.save();
    ctx.setTransform(window.devicePixelRatio || 1, 0, 0, window.devicePixelRatio || 1, 0, 0);
    ctx.clearRect(0, 0, cssW, cssH);
    ctx.fillStyle = '#808080';
    ctx.fillRect(0, 0, cssW, cssH);

    if (phase === 'fixation') {
      drawFixationPhase(ctx, cssW, cssH);
    } else {
      ctx.imageSmoothingEnabled = true;
      ctx.drawImage(trialCanvas, 0, 0, cssW, cssH);
    }
    ctx.restore();
  }

  async function flashTarget(trialCanvas, spec) {
    const startPerf = performance.now();
    const flashUntil = startPerf + state.settings.flashMs;
    return new Promise((resolve) => {
      const step = (frameTs) => {
        const ctx = dom.ctx;
        const cssW = state.stageWidth;
        const cssH = state.stageHeight;
        ctx.save();
        ctx.setTransform(window.devicePixelRatio || 1, 0, 0, window.devicePixelRatio || 1, 0, 0);
        ctx.clearRect(0, 0, cssW, cssH);
        ctx.drawImage(trialCanvas, 0, 0, cssW, cssH);
        drawFlashingRing(ctx, spec, frameTs - startPerf);
        ctx.restore();
        if (frameTs < flashUntil) {
          requestAnimationFrame(step);
        } else {
          resolve(frameTs);
        }
      };
      requestAnimationFrame(step);
    });
  }

  function drawFlashingRing(ctx, spec, elapsedMs) {
    const visible = Math.floor(elapsedMs / 90) % 2 === 0;
    if (!visible) {
      return;
    }
    const x = spec.targetX * state.displayScaleX;
    const y = spec.targetY * state.displayScaleY;
    const radius = Math.max(10, spec.targetRadiusPx * state.displayScaleX * 1.3);
    ctx.save();
    ctx.strokeStyle = 'rgba(255, 239, 191, 0.96)';
    ctx.lineWidth = Math.max(3, Math.round(radius * 0.12));
    ctx.shadowColor = 'rgba(242, 184, 75, 0.28)';
    ctx.shadowBlur = 14;
    ctx.beginPath();
    ctx.arc(x, y, radius, 0, Math.PI * 2);
    ctx.stroke();
    ctx.restore();
  }

  function drawIdleScreen() {
    resizeStageCanvas();
    const ctx = dom.ctx;
    const cssW = state.stageWidth;
    const cssH = state.stageHeight;
    ctx.save();
    ctx.setTransform(window.devicePixelRatio || 1, 0, 0, window.devicePixelRatio || 1, 0, 0);
    ctx.clearRect(0, 0, cssW, cssH);
    ctx.fillStyle = 'rgb(128,128,128)';
    ctx.fillRect(0, 0, cssW, cssH);

    ctx.fillStyle = 'rgba(0, 0, 0, 0.28)';
    ctx.font = '600 22px Georgia, serif';
    ctx.textAlign = 'center';
    ctx.fillText('等待准备实验', cssW / 2, cssH / 2 - 18);
    ctx.font = '16px "Trebuchet MS", sans-serif';
    ctx.fillStyle = 'rgba(0, 0, 0, 0.5)';
    ctx.fillText('准备完成后会自动进入固定中央注视点流程', cssW / 2, cssH / 2 + 18);
    ctx.restore();
  }

  function buildTrialSpecs(settings, sessionSeed, renderWidth, renderHeight) {
    const rng = mulberry32(sessionSeed);
    const specs = [];
    const centerX = renderWidth / 2;
    const centerY = renderHeight / 2;
    const minDim = Math.min(renderWidth, renderHeight);
    const bgRadiusFromHeight = (renderHeight * settings.noiseDiskHeightRatio) / 2;
    const bgRadiusPx = Math.min(bgRadiusFromHeight, minDim * 0.495);
    const pxPerDeg = bgRadiusPx / (settings.bgWidthDeg / 2);
    const targetSigmaPx = Math.max(1, (settings.gaborWidthDeg * pxPerDeg) / 2.355);
    const targetRadiusPx = Math.max(3, targetSigmaPx * 2.4);
    const targetMinEccPx = settings.eccMinDeg * pxPerDeg;
    const targetMaxEccPx = Math.min(settings.eccMaxDeg * pxPerDeg, bgRadiusPx - targetRadiusPx * 1.2);

    for (let trialIndex = 0; trialIndex < settings.trialCount; trialIndex += 1) {
      const hueDeg = rng() * 360;
      const saturation = 0.56 + rng() * 0.08;
      const radiusSample = Math.sqrt(rng());
      const targetEccPx = targetMinEccPx + radiusSample * (targetMaxEccPx - targetMinEccPx);
      const targetOriDeg = rng() * 360;
      const targetX = centerX + targetEccPx * Math.cos(degToRad(targetOriDeg));
      const targetY = centerY + targetEccPx * Math.sin(degToRad(targetOriDeg));
      const fixationHoldMs = lerp(settings.fixMinMs, settings.fixMaxMs, rng());
      const carrierCyclesPerPx = settings.gaborSF / pxPerDeg;

      specs.push({
        trialIndex: trialIndex + 1,
        seed: Math.floor(rng() * 1e9),
        hueDeg,
        saturation,
        bgContrast: settings.bgContrast,
        targetContrast: settings.targetContrast,
        targetX,
        targetY,
        targetOriDeg: settings.gaborOrientDeg,
        targetRadiusPx,
        targetSigmaPx,
        bgRadiusPx,
        fixationHoldMs,
        responseTimeoutMs: settings.responseTimeoutMs,
        carrierCyclesPerPx,
      });
    }
    return specs;
  }

  function buildTrialCanvas(spec, renderWidth, renderHeight) {
    const canvas = document.createElement('canvas');
    canvas.width = renderWidth;
    canvas.height = renderHeight;
    const ctx = canvas.getContext('2d', { willReadFrequently: true });
    const imageData = ctx.createImageData(renderWidth, renderHeight);
    const data = imageData.data;
    const rng = mulberry32(spec.seed);
    const noise = createPinkNoise(renderWidth, renderHeight, rng, spec.bgContrast);
    const rgbScale = hsvValueToRgbScale(spec.hueDeg / 360, spec.saturation);
    const cx = renderWidth / 2;
    const cy = renderHeight / 2;
    const bgRadiusSq = spec.bgRadiusPx * spec.bgRadiusPx;
    const patchRadius = Math.ceil(spec.targetSigmaPx * 4.2);
    const patchX0 = Math.max(0, Math.floor(spec.targetX - patchRadius));
    const patchX1 = Math.min(renderWidth - 1, Math.ceil(spec.targetX + patchRadius));
    const patchY0 = Math.max(0, Math.floor(spec.targetY - patchRadius));
    const patchY1 = Math.min(renderHeight - 1, Math.ceil(spec.targetY + patchRadius));
    const theta = degToRad(spec.targetOriDeg);
    const cosT = Math.cos(theta);
    const sinT = Math.sin(theta);
    const targetSigmaSq2 = 2 * spec.targetSigmaPx * spec.targetSigmaPx;

    for (let y = 0; y < renderHeight; y += 1) {
      const rowOffset = y * renderWidth;
      const dy = y - cy;
      for (let x = 0; x < renderWidth; x += 1) {
        const idx = rowOffset + x;
        const dx = x - cx;
        let value = noise[idx];
        let isNoiseRegion = true;

        if ((dx * dx) + (dy * dy) > bgRadiusSq) {
          value = 0.5;
          isNoiseRegion = false;
        }

        if (x >= patchX0 && x <= patchX1 && y >= patchY0 && y <= patchY1) {
          const tx = x - spec.targetX;
          const ty = y - spec.targetY;
          const xr = tx * cosT + ty * sinT;
          const env = Math.exp(-((tx * tx) + (ty * ty)) / targetSigmaSq2);
          const carrier = Math.cos(2 * Math.PI * spec.carrierCyclesPerPx * xr);
          value = clamp01(value + spec.targetContrast * env * carrier);
        }

        const base = idx * 4;
        if (isNoiseRegion) {
          data[base] = clamp255(value * rgbScale[0] * 255);
          data[base + 1] = clamp255(value * rgbScale[1] * 255);
          data[base + 2] = clamp255(value * rgbScale[2] * 255);
        } else {
          const gray = clamp255(value * 255);
          data[base] = gray;
          data[base + 1] = gray;
          data[base + 2] = gray;
        }
        data[base + 3] = 255;
      }
    }

    ctx.putImageData(imageData, 0, 0);
    return canvas;
  }

  function createPinkNoise(width, height, rng, bgContrast) {
    const weights = [1.15, 0.9, 0.7, 0.52, 0.38, 0.26];
    const scales = [180, 120, 76, 46, 28, 14];
    const layers = weights.map((weight, index) => {
      const gridW = Math.max(4, Math.round(width / scales[index]));
      const gridH = Math.max(4, Math.round(height / scales[index]));
      const values = new Float32Array(gridW * gridH);
      for (let i = 0; i < values.length; i += 1) {
        values[i] = rng();
      }
      return { gridW, gridH, values, weight };
    });

    const output = new Float32Array(width * height);
    let weightSum = 0;
    for (const layer of layers) {
      weightSum += layer.weight;
    }

    for (let y = 0; y < height; y += 1) {
      const yn = height <= 1 ? 0 : y / (height - 1);
      for (let x = 0; x < width; x += 1) {
        const xn = width <= 1 ? 0 : x / (width - 1);
        let sum = 0;
        for (const layer of layers) {
          sum += sampleBilinearLayer(layer.values, layer.gridW, layer.gridH, xn, yn) * layer.weight;
        }
        const normalized = sum / weightSum;
        const contrastBoost = 0.55 + bgContrast * 1.45;
        output[(y * width) + x] = clamp01(0.5 + (normalized - 0.5) * contrastBoost);
      }
    }

    return output;
  }

  function sampleBilinearLayer(values, gridW, gridH, xn, yn) {
    const gx = xn * (gridW - 1);
    const gy = yn * (gridH - 1);
    const x0 = Math.floor(gx);
    const y0 = Math.floor(gy);
    const x1 = Math.min(x0 + 1, gridW - 1);
    const y1 = Math.min(y0 + 1, gridH - 1);
    const tx = gx - x0;
    const ty = gy - y0;
    const v00 = values[(y0 * gridW) + x0];
    const v10 = values[(y0 * gridW) + x1];
    const v01 = values[(y1 * gridW) + x0];
    const v11 = values[(y1 * gridW) + x1];
    const a = v00 + (v10 - v00) * tx;
    const b = v01 + (v11 - v01) * tx;
    return a + (b - a) * ty;
  }

  function hsvValueToRgbScale(hue, saturation) {
    const h = ((hue % 1) + 1) % 1;
    const sector = Math.floor(h * 6);
    const f = h * 6 - sector;
    const p = 1 - saturation;
    const q = 1 - (saturation * f);
    const t = 1 - (saturation * (1 - f));
    switch (sector % 6) {
      case 0:
        return [1, t, p];
      case 1:
        return [q, 1, p];
      case 2:
        return [p, 1, t];
      case 3:
        return [p, q, 1];
      case 4:
        return [t, p, 1];
      default:
        return [1, p, q];
    }
  }

  function waitForResponse(timeoutMs) {
    return new Promise((resolve) => {
      const cleanup = () => {
        if (state.activeResponseTimeoutId !== null) {
          clearTimeout(state.activeResponseTimeoutId);
          state.activeResponseTimeoutId = null;
        }
        state.activeKeyResolver = null;
      };

      state.activeKeyResolver = (payload) => {
        cleanup();
        resolve(payload);
      };

      if (Number.isFinite(timeoutMs) && timeoutMs >= 0) {
        state.activeResponseTimeoutId = window.setTimeout(() => {
          cleanup();
          resolve({
            responded: false,
            timedOut: true,
            perfNow: performance.now(),
            rawTimeStamp: null,
            keyCode: '',
            keyLabel: '',
            eventType: 'timeout',
            timestamp: makeStamp(performance.now()),
          });
        }, timeoutMs);
      }
    });
  }
  
  function handleKeyDown(event) {
      // 忽略自动重复的按键
      if (event.repeat) return;

      /* ---------- 休息态：由 takeBreak 内置的监听器处理 ---------- */
      if (state.breakActive) return;

      /* ===================== 实验正在运行 ====================== */
      if (state.running) {
          // Escape 中止实验
          if (event.code === 'Escape') {
              event.preventDefault();
              abortExperiment('用户按下 Escape 中止。');
              return;
          }

          // 响应按键（主动响应注视点/刺激）
          if (state.activeKeyResolver && event.code === state.settings.responseKey) {
              event.preventDefault();
              const perfNow = performance.now();
              const payload = {
                  responded: true,
                  timedOut: false,
                  perfNow,
                  rawTimeStamp: event.timeStamp,
                  keyCode: event.code,
                  keyLabel: event.key,
                  eventType: event.type,
                  timestamp: makeStamp(perfNow),
              };
              state.activeKeyResolver(payload);
          }
          return;
      }

      /* ===================== 实验未运行 ====================== */
      const overlayVisible = dom.overlay && dom.overlay.style.display !== 'none';

      // 回车键：在向导中推进流程，在主界面启动实验
      if (event.code === 'Enter' && !event.ctrlKey && !event.metaKey) {
          event.preventDefault();

          if (overlayVisible) {
              // 当前处于覆盖层：根据当前显示的卡片模拟点击按钮
              if (dom.infoSection && dom.infoSection.style.display !== 'none') {
                  // 信息卡片 → 进入指导语卡片
                  dom.overlayNextBtn && dom.overlayNextBtn.click();
              } else if (dom.instrSection && dom.instrSection.style.display !== 'none') {
                  // 指导语卡片 → 正式开始实验（含全屏、准备、自动启动）
                  dom.instrStartBtn && dom.instrStartBtn.click();
              }
              // 结束卡片（finishSection）显示时不响应回车，避免误触
          } else if (dom.appShell && dom.appShell.style.display !== 'none') {
              // 主界面已显示 → 允许回车快捷启动实验（原有行为）
              startExperiment().catch(handleFatalError);
          }
          return;
      }

      // Escape 键：在覆盖层中可作为安全退出（也可不处理，由具体卡片决定）
      if (event.code === 'Escape') {
          event.preventDefault();
          if (overlayVisible) {
              // 在信息/指导语卡片按 Escape：不做任何破坏性操作，仅阻止默认行为
              // （若需要可增加关闭覆盖层动作，此处保持温和）
          } else if (dom.appShell && dom.appShell.style.display !== 'none') {
              abortExperiment('用户按下 Escape 中止。');
          }
      }
  }
  function abortExperiment(message) {
    state.aborted = true;
    state.running = false;
    cancelServerShutdown();
    if (state.activeResponseTimeoutId !== null) {
      clearTimeout(state.activeResponseTimeoutId);
      state.activeResponseTimeoutId = null;
    }
    state.activeKeyResolver = null;
    document.body.classList.remove('running');
    exitFullscreenSafe();
    drawIdleScreen();
    setStatus(message || '实验已中止。');
    updateButtons();
  }

  function buildSummaryText() {
    const total = state.results.length;
    const responded = state.results.filter((r) => r.responded).length;
    const timedOut = state.results.filter((r) => r.timeout).length;
    const flashed = state.results.filter((r) => r.flashed).length;
    const rts = state.results.map((r) => r.responseRTMs).filter((v) => typeof v === 'number' && Number.isFinite(v));
    const meanRt = rts.length ? rts.reduce((a, b) => a + b, 0) / rts.length : NaN;
    return [
      `Trials completed: ${total}`,
      `Responses: ${responded}`,
      `Timeouts: ${timedOut}`,
      `Flashes shown: ${flashed}`,
      `Mean RT (ms): ${Number.isFinite(meanRt) ? meanRt.toFixed(1) : 'n/a'}`,
      `Time zone: ${state.settings ? state.settings.timeZone : ''}`,
      `Time origin (ms): ${performance.timeOrigin.toFixed(3)}`,
    ].join('\n');
  }

  function updateLiveSummary(lastTrial) {
    const summary = buildSummaryText();
    const extra = [
      '',
      `Last trial: ${lastTrial.trialIndex}`,
      `Stimulus: ${lastTrial.stimulus ? lastTrial.stimulus.localCompact : 'n/a'}`,
      `Response: ${lastTrial.responded ? `${lastTrial.response.localCompact} | ${lastTrial.responseRTMs.toFixed(1)} ms` : 'timeout'}`,
    ].join('\n');
    setSummary(`${summary}${extra}`);
  }

  function appendTrialLog(result, message) {
    const lines = [
      `Trial ${result.trialIndex}: ${message}`,
      `  fixationPress: ${result.fixationPress ? result.fixationPress.iso : 'n/a'}`,
      `  fixationSolid: ${result.fixationSolid ? result.fixationSolid.iso : 'n/a'}`,
      `  stimulus: ${result.stimulus ? result.stimulus.iso : 'n/a'}`,
    ];
    if (result.responded) {
      lines.push(`  response: ${result.response.iso} (${result.responseRTMs.toFixed(1)} ms)`);
    } else {
      lines.push('  response: timeout');
    }
    appendStatus(lines.join('\n'));
  }

  function downloadResults(kind) {
    if (!state.results.length) {
      return;
    }

    const baseName = buildBaseName();
    let blob;
    let fileName;
    if (kind === 'json') {
      blob = new Blob([JSON.stringify(buildJsonPayload(), null, 2)], { type: 'application/json;charset=utf-8' });
      fileName = `${baseName}.json`;
    } else if (kind === 'csv') {
      const csv = '\ufeff' + buildCsvPayload();
      blob = new Blob([csv], { type: 'text/csv;charset=utf-8' });
      fileName = `${baseName}.csv`;
    } else {
      blob = new Blob([buildTxtPayload()], { type: 'text/plain;charset=utf-8' });
      fileName = `${baseName}.txt`;
    }
    downloadBlob(blob, fileName);
  }

  function buildJsonPayload() {
    return {
      metadata: {
        createdAt: makeStamp(performance.now()),
        sessionStamp: state.sessionStamp,
        settings: state.settings,
        sessionSeed: state.sessionSeed,
        timeOriginMs: performance.timeOrigin,
        browser: {
          userAgent: navigator.userAgent || '',
          language: navigator.language || '',
          platform: navigator.platform || '',
          devicePixelRatio: window.devicePixelRatio || 1,
        },
      },
      trials: state.results,
    };
  }

  function buildCsvPayload() {
    const rows = state.results.map((trial) => ({
      trialIndex: trial.trialIndex,
      seed: trial.seed,
      hueDeg: trial.hueDeg,
      saturation: trial.saturation,
      bgContrast: trial.bgContrast,
      targetContrast: trial.targetContrast,
      fixationHoldMs: trial.fixationHoldMs,
      responseTimeoutMs: trial.responseTimeoutMs,
      renderWidth: trial.renderWidth,
      renderHeight: trial.renderHeight,
      stageWidth: trial.stageWidth,
      stageHeight: trial.stageHeight,
      targetXpx: trial.targetXpx,
      targetYpx: trial.targetYpx,
      targetXnorm: trial.targetXnorm,
      targetYnorm: trial.targetYnorm,
      targetRadiusPx: trial.targetRadiusPx,
      targetSigmaPx: trial.targetSigmaPx,
      targetOrientationDeg: trial.targetOrientationDeg,
      carrierCyclesPerPx: trial.carrierCyclesPerPx,
      bgRadiusPx: trial.bgRadiusPx,
      fixationPress_epochMs: trial.fixationPress ? trial.fixationPress.epochMs : '',
      fixationPress_iso: trial.fixationPress ? trial.fixationPress.iso : '',
      fixationPress_utc: trial.fixationPress ? trial.fixationPress.utcString : '',
      fixationPress_local: trial.fixationPress ? trial.fixationPress.localCompact : '',
      fixationPress_perfMs: trial.fixationPress ? trial.fixationPress.perfMs : '',
      fixationSolid_epochMs: trial.fixationSolid ? trial.fixationSolid.epochMs : '',
      fixationSolid_iso: trial.fixationSolid ? trial.fixationSolid.iso : '',
      fixationSolid_utc: trial.fixationSolid ? trial.fixationSolid.utcString : '',
      fixationSolid_local: trial.fixationSolid ? trial.fixationSolid.localCompact : '',
      fixationSolid_perfMs: trial.fixationSolid ? trial.fixationSolid.perfMs : '',
      stimulus_epochMs: trial.stimulus ? trial.stimulus.epochMs : '',
      stimulus_iso: trial.stimulus ? trial.stimulus.iso : '',
      stimulus_utc: trial.stimulus ? trial.stimulus.utcString : '',
      stimulus_local: trial.stimulus ? trial.stimulus.localCompact : '',
      stimulus_perfMs: trial.stimulus ? trial.stimulus.perfMs : '',
      response_epochMs: trial.response ? trial.response.epochMs : '',
      response_iso: trial.response ? trial.response.iso : '',
      response_utc: trial.response ? trial.response.utcString : '',
      response_local: trial.response ? trial.response.localCompact : '',
      response_perfMs: trial.response ? trial.response.perfMs : '',
      responseRTMs: trial.responseRTMs ?? '',
      responseKey: trial.responseKey || '',
      responseKeyLabel: trial.responseKeyLabel || '',
      responseTimedOut: trial.timeout,
      questJudge: trial.questJudge ?? '',
      questRtCriterionMs: trial.questRtCriterionMs ?? '',
      questTargetP: trial.questTargetP ?? '',
      questGrain: trial.questGrain ?? '',
      questBoostLog10: trial.questBoostLog10 ?? '',
      questIntensityLog10: trial.questIntensityLog10 ?? '',
      backgroundFlash_epochMs: trial.backgroundFlash ? trial.backgroundFlash.epochMs : '',
      backgroundFlash_iso: trial.backgroundFlash ? trial.backgroundFlash.iso : '',
      flash_epochMs: trial.flash ? trial.flash.epochMs : '',
      flash_iso: trial.flash ? trial.flash.iso : '',
      flash_utc: trial.flash ? trial.flash.utcString : '',
      flash_local: trial.flash ? trial.flash.localCompact : '',
      flash_perfMs: trial.flash ? trial.flash.perfMs : '',
      trialEnd_epochMs: trial.trialEnd ? trial.trialEnd.epochMs : '',
      trialEnd_iso: trial.trialEnd ? trial.trialEnd.iso : '',
      trialEnd_utc: trial.trialEnd ? trial.trialEnd.utcString : '',
      trialEnd_local: trial.trialEnd ? trial.trialEnd.localCompact : '',
      trialEnd_perfMs: trial.trialEnd ? trial.trialEnd.perfMs : '',
      timeZone: state.settings.timeZone,
      locale: state.settings.locale,
      timeOriginMs: performance.timeOrigin,
      sessionStart_epochMs: state.sessionStamp ? state.sessionStamp.epochMs : '',
      sessionStart_iso: state.sessionStamp ? state.sessionStamp.iso : '',
      sessionStart_utc: state.sessionStamp ? state.sessionStamp.utcString : '',
      sessionStart_local: state.sessionStamp ? state.sessionStamp.localCompact : '',
      sessionStart_perfMs: state.sessionStamp ? state.sessionStamp.perfMs : '',
      subjectId: state.settings.subjectId,
      subjectName: state.settings.subjectName,
      subjectGender: state.settings.subjectGender,
      organization: state.settings.organization,
      userAgent: state.settings.userAgent,
      screenWidth: state.settings.screenWidth,
      screenHeight: state.settings.screenHeight,
      pageWidth: state.settings.pageWidth,
      pageHeight: state.settings.pageHeight,
      devicePixelRatio: state.settings.devicePixelRatio,
    }));
    return toCsv(rows);
  }

  function buildTxtPayload() {
    const lines = [];
    lines.push('FreeView Web Experiment');
    lines.push(`Session seed: ${state.sessionSeed}`);
    lines.push(`Session start: ${state.sessionStamp ? state.sessionStamp.iso : ''}`);
    lines.push(`Participant: ${state.settings.subjectId}, ${state.settings.subjectName}, ${state.settings.subjectGender}, ${state.settings.organization}`);
    lines.push(`Time origin: ${performance.timeOrigin.toFixed(3)} ms`);
    lines.push('');
    for (const trial of state.results) {
      lines.push(`Trial ${trial.trialIndex}`);
      lines.push(`  seed: ${trial.seed}`);
      lines.push(`  fixationPress: ${trial.fixationPress ? trial.fixationPress.iso : 'n/a'}`);
      lines.push(`  fixationSolid: ${trial.fixationSolid ? trial.fixationSolid.iso : 'n/a'}`);
      lines.push(`  stimulus: ${trial.stimulus ? trial.stimulus.iso : 'n/a'}`);
      lines.push(`  response: ${trial.responded ? `${trial.response.iso} (${trial.responseRTMs.toFixed(1)} ms)` : 'timeout'}`);
      lines.push(`  questJudge: ${trial.questJudge} (criterion ${trial.questRtCriterionMs} ms)`);
      lines.push(`  flash: ${trial.flash ? trial.flash.iso : 'n/a'}`);
      lines.push(`  target: x=${trial.targetXpx.toFixed(2)}, y=${trial.targetYpx.toFixed(2)}`);
      lines.push('');
    }
    return `${lines.join('\n')}\n`;
  }

  function toCsv(rows) {
    if (!rows.length) {
      return '';
    }
    const headers = Object.keys(rows[0]);
    const lines = [headers.join(',')];
    for (const row of rows) {
      lines.push(headers.map((header) => csvEscape(row[header])).join(','));
    }
    return `${lines.join('\n')}\n`;
  }

  function csvEscape(value) {
    if (value === null || value === undefined || Number.isNaN(value)) {
      return '';
    }
    const text = String(value);
    if (/[",\n\r]/.test(text)) {
      return `"${text.replace(/"/g, '""')}"`;
    }
    return text;
  }

  function buildBaseName() {
    const subject = sanitizeFilePart(state.settings ? state.settings.subjectId : 'S001');
    const session = sanitizeFilePart(state.settings ? state.settings.sessionId : '1');
    const stamp = formatCompactFileStamp(new Date());
    return `WebFreeView_${subject}_Ses${session}_${stamp}`;
  }

  function downloadBlob(blob, fileName) {
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = fileName;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
  }

  function updateButtons() {
    dom.startBtn.disabled = !state.prepared || state.running;
  }

  function scheduleAutoPrepare(delayMs) {
    if (state.running) {
      return;
    }
    state.prepared = false;
    updateButtons();
    if (state.autoPrepareTimerId !== null) {
      window.clearTimeout(state.autoPrepareTimerId);
      state.autoPrepareTimerId = null;
    }
    state.autoPrepareTimerId = window.setTimeout(() => {
      state.autoPrepareTimerId = null;
      setStatus('页面参数变更，正在自动准备...');
      prepareExperiment().catch(handleFatalError);
    }, Math.max(0, delayMs));
  }

  function scheduleServerShutdown(delayMs) {
    cancelServerShutdown();
    state.serverShutdownTimerId = window.setTimeout(() => {
      state.serverShutdownTimerId = null;
      requestServerShutdown();
    }, Math.max(0, delayMs));
  }

  function cancelServerShutdown() {
    if (state.serverShutdownTimerId !== null) {
      window.clearTimeout(state.serverShutdownTimerId);
      state.serverShutdownTimerId = null;
    }
  }

  function requestServerShutdown() {
    if (window.location.protocol !== 'http:' && window.location.protocol !== 'https:') {
      return;
    }
    const shutdownUrl = new URL('__shutdown__', window.location.href).toString();
    fetch(shutdownUrl, {
      method: 'GET',
      cache: 'no-store',
      keepalive: true,
      credentials: 'same-origin',
    }).catch(() => {});
  }

  function initializeQuestState(settings, trialCount) {
    const learnTNum = trialCount;
    const base = createQuest(
      Math.log10(Math.max(0.001, settings.targetContrast)),
      1,
      settings.learnP,
      3,
      0.01,
      2 * 2 / (settings.bgWidthDeg * settings.bgWidthDeg),
      0.01,
      3
    );
    updateQuestPosterior(base, base.tGuess, 1);
    updateQuestPosterior(base, Math.log10(0.005), 0);
    return {
      q: base,
      learnP: settings.learnP,
      testP: settings.testP,
      connectTNum: settings.connectTNum,
      learnTNum,
      minLog10: Math.log10(0.001),
      maxLog10: Math.log10(1.0),
    };
  }

  function createQuest(tGuess, tGuessSd, pThreshold, beta, delta, gamma, grain, range) {
    const xMin = tGuess - range / 2;
    const xMax = tGuess + range / 2;
    const count = Math.max(2, Math.floor((xMax - xMin) / grain) + 1);
    const grid = new Float64Array(count);
    const posterior = new Float64Array(count);
    let sum = 0;
    for (let i = 0; i < count; i += 1) {
      const t = xMin + i * grain;
      grid[i] = t;
      const z = (t - tGuess) / Math.max(1e-6, tGuessSd);
      const p = Math.exp(-0.5 * z * z);
      posterior[i] = p;
      sum += p;
    }
    for (let i = 0; i < count; i += 1) {
      posterior[i] /= sum;
    }
    const k = Math.max(1e-6, (pThreshold - gamma) / Math.max(1e-6, 1 - gamma - delta));
    const logisticShift = Math.log(k / Math.max(1e-6, 1 - k)) / Math.max(1e-6, beta);
    return {
      grid,
      posterior,
      tGuess,
      pThreshold,
      beta,
      delta,
      gamma,
      grain,
      range,
      logisticShift,
    };
  }

  function computeQuestTrialPlan(questState, history, trialNumber, trialCount) {
    const q = questState.q;
    if (trialNumber > questState.learnTNum && trialNumber <= questState.learnTNum + questState.connectTNum) {
      const step = (questState.learnP - questState.testP) / Math.max(1, questState.connectTNum);
      q.pThreshold = q.pThreshold - step;
      const k = Math.max(1e-6, (q.pThreshold - q.gamma) / Math.max(1e-6, 1 - q.gamma - q.delta));
      q.logisticShift = Math.log(k / Math.max(1e-6, 1 - k)) / Math.max(1e-6, q.beta);
    }

    const recent = history.slice(Math.max(0, history.length - 12));
    const acc = recent.length ? recent.reduce((a, b) => a + b, 0) / recent.length : 0.5;
    let errStreak = 0;
    for (let i = recent.length - 1; i >= 0; i -= 1) {
      if (recent[i] === 0) {
        errStreak += 1;
      } else {
        break;
      }
    }

    let grain;
    let boostLog10;
    if (acc < 0.55 || errStreak >= 3) {
      grain = 0.05;
      boostLog10 = Math.log10(1.35);
    } else if (acc > 0.85 && errStreak === 0) {
      grain = 0.01;
      boostLog10 = Math.log10(0.93);
    } else {
      grain = 0.02;
      boostLog10 = 0;
    }
    q.grain = grain;

    const qMid = questQuantile(q, 0.5);
    const raw = clampFloat(qMid + boostLog10, questState.minLog10, questState.maxLog10, qMid);
    const snapped = Math.round(raw / grain) * grain;
    const tTest = clampFloat(snapped, questState.minLog10, questState.maxLog10, raw);
    return {
      grain,
      boostLog10,
      tTest,
      pThreshold: q.pThreshold,
      targetContrast: Math.pow(10, tTest),
    };
  }

  function getQuestJudgeFromTrial(trialResult, rtCriterionMs) {
    if (!trialResult.responded || trialResult.timeout || typeof trialResult.responseRTMs !== 'number') {
      return 0;
    }
    return trialResult.responseRTMs <= rtCriterionMs ? 1 : 0;
  }

  function updateQuestPosterior(q, intensityLog10, responseCorrect) {
    if (!q || !q.grid || !q.posterior || typeof q.posterior.length !== 'number') {
      return;
    }
    const eps = 1e-9;
    let sum = 0;
    for (let i = 0; i < q.grid.length; i += 1) {
      const p = psychometricProbability(q, intensityLog10, q.grid[i]);
      const likelihood = responseCorrect ? p : (1 - p);
      q.posterior[i] *= Math.max(eps, likelihood);
      sum += q.posterior[i];
    }
    if (sum <= 0) {
      const uni = 1 / q.posterior.length;
      for (let i = 0; i < q.posterior.length; i += 1) {
        q.posterior[i] = uni;
      }
      return;
    }
    for (let i = 0; i < q.posterior.length; i += 1) {
      q.posterior[i] /= sum;
    }
  }

  function psychometricProbability(q, intensity, threshold) {
    const z = q.beta * (intensity - (threshold + q.logisticShift));
    const s = 1 / (1 + Math.exp(-z));
    const p = q.gamma + (1 - q.gamma - q.delta) * s;
    return Math.max(1e-6, Math.min(1 - 1e-6, p));
  }

  function questQuantile(q, percentile) {
    let cumulative = 0;
    const target = Math.max(0, Math.min(1, percentile));
    for (let i = 0; i < q.grid.length; i += 1) {
      cumulative += q.posterior[i];
      if (cumulative >= target) {
        return q.grid[i];
      }
    }
    return q.grid[q.grid.length - 1];
  }

  function setStatus(text) {
    dom.statusText.textContent = text;
  }

  function appendStatus(text) {
    const existing = dom.statusText.textContent || '';
    const lines = existing.split('\n').filter(Boolean);
    const newLines = text.split('\n').filter(Boolean);
    const merged = lines.concat(newLines).slice(-18);
    dom.statusText.textContent = merged.join('\n');
  }

  function setSummary(text) {
    dom.summaryText.textContent = text;
  }

  function setProgress(fraction) {
    const pct = Math.max(0, Math.min(100, fraction * 100));
    dom.prepProgressFill.style.width = `${pct}%`;
  }

  function handleFatalError(error) {
    console.error(error);
    state.running = false;
    state.aborted = true;
    document.body.classList.remove('running');
    setStatus(`错误：${error && error.message ? error.message : String(error)}`);
    setSummary('发生错误，未生成完整结果。');
    drawIdleScreen();
    updateButtons();
  }

  function resizeStageCanvas() {
    const rect = dom.stageFrame.getBoundingClientRect();
    state.stageWidth = Math.max(1, Math.floor(rect.width));
    state.stageHeight = Math.max(1, Math.floor(rect.height));
    const dpr = window.devicePixelRatio || 1;
    dom.canvas.width = Math.max(1, Math.round(state.stageWidth * dpr));
    dom.canvas.height = Math.max(1, Math.round(state.stageHeight * dpr));
    dom.canvas.style.width = `${state.stageWidth}px`;
    dom.canvas.style.height = `${state.stageHeight}px`;
  }

  function rebuildGeometryForCurrentStage() {
    resizeStageCanvas();
    const available = getStageCssSize();
    const dpr = window.devicePixelRatio || 1;
    state.renderWidth = Math.max(1, Math.round(available.width * dpr));
    state.renderHeight = Math.max(1, Math.round(available.height * dpr));
    state.displayScaleX = available.width / state.renderWidth;
    state.displayScaleY = available.height / state.renderHeight;
    state.trialSpecs = buildTrialSpecs(state.settings, state.sessionSeed, state.renderWidth, state.renderHeight);
  }

  function initSubjectCounter() {
    const stored = window.localStorage.getItem('webfreeview_next_subject');
    const parsed = parseInt(stored || '', 10);
    const nextId = Number.isFinite(parsed) && parsed >= 1 ? parsed : 1;
    dom.subjectId.value = String(nextId);
  }

  function incrementSubjectCounter() {
    const current = parseInt(dom.subjectId.value, 10);
    if (!Number.isFinite(current) || current < 1) {
      return;
    }
    const next = current + 1;
    window.localStorage.setItem('webfreeview_next_subject', String(next));
    dom.subjectId.value = String(next);
  }

  function autoDownloadAllResults() {
    if (state.autoDownloadDone || !state.results.length) {
      return;
    }
    state.autoDownloadDone = true;
    // 默认仅自动下载 JSON；CSV/TXT 保留手动下载按钮
    downloadResults('json');
  }

  function copyJsonToClipboard() {
    try {
      const payload = buildJsonPayload();
      const text = JSON.stringify(payload, null, 2);
      if (navigator.clipboard && navigator.clipboard.writeText) {
        navigator.clipboard.writeText(text).then(() => {
          setStatus('JSON 已复制到剪贴板。');
        }).catch(() => {
          setStatus('复制到剪贴板失败。');
        });
      } else {
        const ta = document.createElement('textarea');
        ta.value = text;
        document.body.appendChild(ta);
        ta.select();
        try {
          document.execCommand('copy');
          setStatus('JSON 已复制到剪贴板。');
        } catch (e) {
          setStatus('复制到剪贴板失败。');
        }
        ta.remove();
      }
    } catch (err) {
      console.error(err);
      setStatus('复制到剪贴板失败。');
    }
  }

  function getStageCssSize() {
    const rect = dom.stageFrame.getBoundingClientRect();
    return { width: Math.max(1, Math.floor(rect.width)), height: Math.max(1, Math.floor(rect.height)) };
  }

  function requestFullscreenSafe() {
    if (document.fullscreenElement) {
      return Promise.resolve(true);
    }
    if (!document.documentElement.requestFullscreen) {
      return Promise.resolve(false);
    }
    return document.documentElement.requestFullscreen().then(() => true).catch(() => false);
  }

  function exitFullscreenSafe() {
    if (!document.fullscreenElement || !document.exitFullscreen) {
      return Promise.resolve(false);
    }
    return document.exitFullscreen().then(() => true).catch(() => false);
  }

  function armAutoFullscreen() {
    let attempted = false;
    const trigger = () => {
      if (attempted) return;
      attempted = true;
      requestFullscreenSafe();
      document.removeEventListener('pointerdown', trigger, true);
      document.removeEventListener('keydown', trigger, true);
    };
    document.addEventListener('pointerdown', trigger, true);
    document.addEventListener('keydown', trigger, true);
  }

  function showFinishOverlay() {
    if (!dom.overlay || !dom.finishSection) {
      return;
    }
    dom.infoSection.style.display = 'none';
    dom.instrSection.style.display = 'none';
    dom.finishSection.style.display = 'block';
    dom.overlay.style.display = 'flex';
  }

  function makeStamp(perfMs) {
    const epochMs = performance.timeOrigin + perfMs;
    const date = new Date(epochMs);
    const localCompact = formatLocalDate(date);
    return {
      perfMs: round3(perfMs),
      epochMs: Math.round(epochMs),
      epochSeconds: round6(epochMs / 1000),
      iso: date.toISOString(),
      utcString: date.toUTCString(),
      localCompact,
      localString: date.toLocaleString(),
      dateOnly: formatDateOnly(date),
      timeOnly: formatTimeOnly(date),
      timezoneOffsetMin: date.getTimezoneOffset(),
      timeZone: Intl.DateTimeFormat().resolvedOptions().timeZone || '',
    };
  }

  function formatLocalDate(date) {
    return `${pad(date.getFullYear(), 4)}-${pad(date.getMonth() + 1, 2)}-${pad(date.getDate(), 2)} ${pad(date.getHours(), 2)}:${pad(date.getMinutes(), 2)}:${pad(date.getSeconds(), 2)}.${pad(date.getMilliseconds(), 3)}`;
  }

  function formatDateOnly(date) {
    return `${pad(date.getFullYear(), 4)}-${pad(date.getMonth() + 1, 2)}-${pad(date.getDate(), 2)}`;
  }

  function formatTimeOnly(date) {
    return `${pad(date.getHours(), 2)}:${pad(date.getMinutes(), 2)}:${pad(date.getSeconds(), 2)}.${pad(date.getMilliseconds(), 3)}`;
  }

  function formatCompactFileStamp(date) {
    return `${pad(date.getFullYear(), 4)}${pad(date.getMonth() + 1, 2)}${pad(date.getDate(), 2)}T${pad(date.getHours(), 2)}${pad(date.getMinutes(), 2)}${pad(date.getSeconds(), 2)}`;
  }

  function sanitizeFilePart(value) {
    return String(value || '')
      .replace(/[^0-9A-Za-z_-]+/g, '_')
      .replace(/^_+|_+$/g, '') || 'NA';
  }

  function createSessionSeed(settings) {
    const base = Date.now() ^ Math.floor(Math.random() * 0x7fffffff);
    const subjectHash = hashString(settings.subjectId || '');
    const sessionHash = hashString(settings.sessionId || '');
    return (base ^ subjectHash ^ (sessionHash << 1)) >>> 0;
  }

  function hashString(text) {
    let hash = 2166136261;
    for (let i = 0; i < text.length; i += 1) {
      hash ^= text.charCodeAt(i);
      hash = Math.imul(hash, 16777619);
    }
    return hash >>> 0;
  }

  function mulberry32(seed) {
    let t = seed >>> 0;
    return function () {
      t += 0x6D2B79F5;
      let r = Math.imul(t ^ (t >>> 15), 1 | t);
      r ^= r + Math.imul(r ^ (r >>> 7), 61 | r);
      return ((r ^ (r >>> 14)) >>> 0) / 4294967296;
    };
  }

  function clampInt(value, min, max, fallback) {
    const num = Number.isFinite(value) ? value : fallback;
    return Math.max(min, Math.min(max, Math.round(num)));
  }

  function clampFloat(value, min, max, fallback) {
    const num = Number.isFinite(value) ? value : fallback;
    return Math.max(min, Math.min(max, num));
  }

  function clamp01(value) {
    return Math.max(0, Math.min(1, value));
  }

  function clamp255(value) {
    return Math.max(0, Math.min(255, Math.round(value)));
  }

  function round3(value) {
    return Math.round(value * 1000) / 1000;
  }

  function round4(value) {
    return Math.round(value * 10000) / 10000;
  }

  function round6(value) {
    return Math.round(value * 1000000) / 1000000;
  }

  function lerp(min, max, fraction) {
    return min + (max - min) * fraction;
  }

  function degToRad(deg) {
    return (deg * Math.PI) / 180;
  }

  function pad(value, length) {
    return String(value).padStart(length, '0');
  }

  function nextAnimationFrame() {
    return new Promise((resolve) => {
      requestAnimationFrame((ts) => resolve(ts));
    });
  }

  function waitMs(ms) {
    return new Promise((resolve) => {
      window.setTimeout(resolve, ms);
    });
  }

  function waitUntilFrame(targetPerfMs, drawFn) {
    return new Promise((resolve) => {
      const step = (frameTs) => {
        if (frameTs >= targetPerfMs) {
          drawFn();
          resolve(frameTs);
        } else {
          requestAnimationFrame(step);
        }
      };
      requestAnimationFrame(step);
    });
  }
})();