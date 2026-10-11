'use strict';
let applet, steps, selected=1, ready=false;
const el=id=>document.getElementById(id);
const script=code=>{if(applet)Jmol.script(applet,code);};
function setEnabled(enabled){for(const id of ['play','pause','restart','product','reset','progress'])el(id).disabled=!enabled;}
function updateFrame(n){
  if(!steps)return;
  n=Math.max(1,Math.min(61,Math.round(Number(n)||1)));
  el('progress').value=n;el('frame').textContent=`${n} / 61`;
  el('phase').textContent=steps[selected-1].phases[n<10?0:n>52?2:1];
}
window.onJmolFrame=function(_id,_index,_file,frame){updateFrame(frame);};
window.onJmolMessage=function(_id,msg){
  if(String(msg).includes('ESOMEPRAZOLE_READY')){el('loading').hidden=true;setEnabled(true);updateFrame(1);}
  else if(/script ERROR|Error reading|FileNotFound/i.test(String(msg))){el('loading').hidden=false;el('loading').textContent='The model could not be loaded. Reload this card to retry.';setEnabled(false);}
};
function loadStep(n){
  if(!ready||!steps)return;
  selected=n;const step=steps[n-1];
  el('stepTitle').textContent=`${n}. ${step.title}`;el('conditions').textContent=step.conditions;
  el('description').textContent=step.description;
  el('scope').textContent=n===9?'The detailed titanium catalyst cycle is not represented.':'';
  for(const b of el('steps').children)b.setAttribute('aria-current',String(Number(b.dataset.step)===n));
  setEnabled(false);el('loading').hidden=false;el('loading').textContent='Loading reaction frames…';
  const url=new URL(step.script,document.baseURI).href;
  script(`animation off; script ${JSON.stringify(url)}; save orientation "initial"; set showHydrogens ${el('hydrogens').checked}; print "ESOMEPRAZOLE_READY";`);
}
async function init(){
  try{
    if(!window.Jmol)throw new Error('JSmol library unavailable');
    const response=await fetch('steps.json');if(!response.ok)throw new Error('Model list unavailable');steps=await response.json();
    const number=Number(new URLSearchParams(location.search).get('step'));
    selected=Number.isInteger(number)&&number>=1&&number<=9?number:1;
    for(const step of steps){const b=document.createElement('button');b.type='button';b.textContent=String(step.num);b.title=step.title;b.setAttribute('aria-label',`${step.num}. ${step.title}`);b.dataset.step=step.num;b.addEventListener('click',()=>loadStep(step.num));el('steps').appendChild(b);}
    const info={width:'100%',height:'100%',use:'HTML5',color:'#FFFFFF',j2sPath:'https://chemapps.stolaf.edu/jmol/jsmol/j2s',serverURL:'https://chemapps.stolaf.edu/jmol/jsmol/php/jsmol.php',allowJavaScript:false,disableInitialConsole:true,disableJ2SLoadMonitor:true,animFrameCallback:'onJmolFrame',messageCallback:'onJmolMessage',readyFunction:function(a){applet=a;ready=true;loadStep(selected);}};
    Jmol.setDocument(0);el('viewer').innerHTML=Jmol.getAppletHtml('esomeprazoleJSmol',info);
    el('play').onclick=()=>script(`animation mode ${el('loop').checked?'LOOP 1 1':'ONCE'}; frame ${Number(el('progress').value)>=61?'1; frame play':'play'};`);
    el('pause').onclick=()=>script('frame pause;');
    el('restart').onclick=()=>script('animation off; frame 1;');
    el('product').onclick=()=>script('animation off; frame 61;');
    el('reset').onclick=()=>script('restore orientation "initial" 0;');
    el('progress').oninput=()=>{const n=Number(el('progress').value);script(`animation off; frame ${n};`);updateFrame(n);};
    el('hydrogens').onchange=()=>script(`set showHydrogens ${el('hydrogens').checked};`);
    el('loop').onchange=()=>script(`animation mode ${el('loop').checked?'LOOP 1 1':'ONCE'};`);
  }catch(error){el('loading').textContent='Unable to initialise JSmol. Reload this card to retry.';console.error(error);}
}
init();
