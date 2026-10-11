'use strict';
let applet,steps,selected=1,ready=false,loadingModel=false,diagramScale=1,diagramSource='',overview=false,extra=false;
const VERSION='paired-20261011b';
const el=id=>document.getElementById(id);
const current=()=>steps[selected-1];
const script=code=>{if(applet)Jmol.script(applet,code);};
function setEnabled(enabled){for(const id of ['play','pause','restart','product','reset','zoomOut','zoomIn','progress'])el(id).disabled=!enabled;for(const b of el('steps').children)b.disabled=!enabled;for(const b of el('stages').children)b.disabled=!enabled;}
function stageFor(n){return current().stages?.find(s=>n>=s.start&&n<=s.end);}
function displayDiagram(src,caption){
  if(src!==diagramSource){diagramSource=src;diagramScale=1;el('diagram').style.width='100%';el('diagram').src=src;el('diagramViewport').scrollTop=0;el('diagramViewport').scrollLeft=0;}
  el('diagram').alt=caption;el('diagramCaption').textContent=caption;el('downloadDiagram').href=src;el('downloadDiagram').download=src.split('/').pop();
}
function updateFrame(n){
  if(!steps)return;n=Math.max(1,Math.min(current().frames,Math.round(Number(n)||1)));
  el('progress').value=n;el('frame').textContent=`${n} / ${current().frames}`;
  const stage=stageFor(n);el('phase').textContent=stage?stage.name:current().phases[n<10?0:n>current().frames-9?2:1];
  for(const b of el('stages').children)b.setAttribute('aria-current',String(stage&&Number(b.dataset.frame)===stage.start));
  if(overview)displayDiagram('diagrams/overview.png','Complete ChemDraw synthesis route');
  else if(extra&&current().extraDiagram)displayDiagram(current().extraDiagram.src,current().extraDiagram.name);
  else displayDiagram(stage?stage.diagram:current().diagram,stage?stage.name:`Step ${selected} — native ChemDraw mechanism`);
}
window.onJmolFrame=function(_id,_index,_file,frame){updateFrame(frame);};
window.onJmolMessage=function(_id,msg){
  if(String(msg).includes(`ESOMEPRAZOLE_READY_${selected}`)){loadingModel=false;el('loading').hidden=true;setEnabled(true);updateFrame(1);}
  else if(/script ERROR|Error reading|FileNotFound/i.test(String(msg))){console.error(msg);if(loadingModel){el('loading').hidden=false;el('loading').textContent='The model could not be loaded. Reload this card to retry.';setEnabled(false);}}
};
function loadStep(n){
  if(!ready||!steps)return;selected=n;const step=current();overview=false;extra=false;
  el('showOverview').setAttribute('aria-pressed','false');el('showExtra').setAttribute('aria-pressed','false');el('showExtra').hidden=!step.extraDiagram;
  el('stepTitle').textContent=`${n}. ${step.title}`;el('conditions').textContent=step.conditions;el('description').textContent=step.description;
  el('progress').max=step.frames;el('stages').replaceChildren();el('stages').hidden=!step.stages;
  for(const stage of step.stages||[]){const b=document.createElement('button');b.type='button';b.textContent=stage.name;b.dataset.frame=stage.start;b.onclick=()=>{overview=false;extra=false;el('showOverview').setAttribute('aria-pressed','false');script(`animation off; frame ${stage.start};`);updateFrame(stage.start);};el('stages').appendChild(b);}
  for(const b of el('steps').children)b.setAttribute('aria-current',String(Number(b.dataset.step)===n));
  updateFrame(1);loadingModel=true;setEnabled(false);el('loading').hidden=false;el('loading').textContent='Loading reaction frames…';
  const url=new URL(step.script,document.baseURI);url.searchParams.set('v',VERSION);
  script(`animation off; script ${JSON.stringify(url.href)}; save orientation "initial"; set showHydrogens ${el('hydrogens').checked}; print "ESOMEPRAZOLE_READY_${n}";`);
}
function zoomDiagram(factor){diagramScale=Math.max(1,Math.min(4,diagramScale*factor));el('diagram').style.width=`${diagramScale*100}%`;}
async function init(){
  try{
    if(!window.Jmol)throw new Error('JSmol library unavailable');
    const response=await fetch(`steps.json?v=${VERSION}`);if(!response.ok)throw new Error('Model list unavailable');steps=await response.json();
    const number=Number(new URLSearchParams(location.search).get('step'));selected=Number.isInteger(number)&&number>=1&&number<=9?number:1;
    for(const step of steps){const b=document.createElement('button');b.type='button';b.textContent=String(step.num);b.title=step.title;b.setAttribute('aria-label',`${step.num}. ${step.title}`);b.dataset.step=step.num;b.onclick=()=>loadStep(step.num);el('steps').appendChild(b);}
    const info={width:'100%',height:'100%',use:'HTML5',color:'#FFFFFF',j2sPath:'https://chemapps.stolaf.edu/jmol/jsmol/j2s',serverURL:'https://chemapps.stolaf.edu/jmol/jsmol/php/jsmol.php',allowJavaScript:false,disableInitialConsole:true,disableJ2SLoadMonitor:true,animFrameCallback:'onJmolFrame',messageCallback:'onJmolMessage',readyFunction:function(a){applet=a;ready=true;loadStep(selected);}};
    Jmol.setDocument(0);el('viewer').innerHTML=Jmol.getAppletHtml('esomeprazoleJSmol',info);
    el('play').onclick=()=>script(`animation mode ${el('loop').checked?'LOOP 1 1':'ONCE'}; frame ${Number(el('progress').value)>=current().frames?'1; frame play':'play'};`);
    el('pause').onclick=()=>script('frame pause;');el('restart').onclick=()=>script('animation off; frame 1;');el('product').onclick=()=>script(`animation off; frame ${current().frames};`);
    el('reset').onclick=()=>script('restore orientation "initial" 0;');el('zoomIn').onclick=()=>script('zoom *1.2;');el('zoomOut').onclick=()=>script('zoom /1.2;');
    el('progress').oninput=()=>{const n=Number(el('progress').value);script(`animation off; frame ${n};`);updateFrame(n);};
    el('hydrogens').onchange=()=>script(`set showHydrogens ${el('hydrogens').checked};`);el('loop').onchange=()=>script(`animation mode ${el('loop').checked?'LOOP 1 1':'ONCE'};`);
    el('showOverview').onclick=()=>{overview=!overview;extra=false;el('showOverview').setAttribute('aria-pressed',String(overview));el('showExtra').setAttribute('aria-pressed','false');updateFrame(el('progress').value);};
    el('showExtra').onclick=()=>{extra=!extra;overview=false;el('showExtra').setAttribute('aria-pressed',String(extra));el('showOverview').setAttribute('aria-pressed','false');updateFrame(el('progress').value);};
    el('diagramIn').onclick=()=>zoomDiagram(1.35);el('diagramOut').onclick=()=>zoomDiagram(1/1.35);el('diagramFit').onclick=()=>{diagramScale=1;el('diagram').style.width='100%';};
    el('expandDiagram').onclick=()=>{script('frame pause;');el('largeDiagram').src=diagramSource;el('dialogTitle').textContent=el('diagramCaption').textContent;el('diagramDialog').showModal();};el('closeDiagram').onclick=()=>el('diagramDialog').close();
    let drag;
    el('diagramViewport').onpointerdown=e=>{if(diagramScale<=1)return;drag={x:e.clientX,y:e.clientY,left:el('diagramViewport').scrollLeft,top:el('diagramViewport').scrollTop};el('diagramViewport').setPointerCapture(e.pointerId);el('diagramViewport').classList.add('dragging');};
    el('diagramViewport').onpointermove=e=>{if(!drag)return;el('diagramViewport').scrollLeft=drag.left+drag.x-e.clientX;el('diagramViewport').scrollTop=drag.top+drag.y-e.clientY;};
    el('diagramViewport').onpointerup=()=>{drag=null;el('diagramViewport').classList.remove('dragging');};
  }catch(error){el('loading').textContent='Unable to initialise JSmol. Reload this card to retry.';console.error(error);}
}
init();
