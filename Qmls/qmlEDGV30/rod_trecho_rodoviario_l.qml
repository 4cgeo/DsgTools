<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'><qgis version="2.6.0-Brighton" minimumScale="1" maximumScale="1" simplifyDrawingHints="0" minLabelScale="0" maxLabelScale="1e+08" simplifyDrawingTol="1" simplifyMaxScale="1" hasScaleBasedVisibilityFlag="0" simplifyLocal="1" scaleBasedLabelVisibilityFlag="0"> 
  <edittypes> 
     <edittype widgetv2type="TextEdit" name="OGC_FID"> 
      <widgetv2config IsMultiline="0" fieldEditable="0" UseHtml="0" labelOnTop="0"/> 
    </edittype> 
    <edittype widgetv2type="TextEdit" name="id"> 
      <widgetv2config IsMultiline="0" fieldEditable="0" UseHtml="0" labelOnTop="0"/> 
    </edittype>
    <edittype widgetv2type="ValueMap" name="geometriaaproximada">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="N�o" value="0"/>
        <value key="Sim" value="1"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="acostamento">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="N�o" value="0"/>
        <value key="Sim" value="1"/>
        <value key="Auto-estrada" value="1"/>
        <value key="Rodovia" value="3"/>
        <value key="Entroncamento" value="4"/>
        <value key="Internacional" value="1"/>
        <value key="Propriedade particular" value="2"/>
        <value key="Federal" value="3"/>
        <value key="Estadual/Distrital" value="4"/>
        <value key="Municipal" value="5"/>
        <value key="Desconhecida" value="95"/>
        <value key="Federal" value="1"/>
        <value key="Estadual/Distrital" value="2"/>
        <value key="Municipal" value="3"/>
        <value key="Concessionada" value="4"/>
        <value key="Privada" value="5"/>
        <value key="N�o aplic�vel" value="6"/>
        <value key="Desconhecida" value="95"/>
        <value key="Revestimento prim�rio (solto)" value="1"/>
        <value key="Pavimentado" value="2"/>
        <value key="Madeira" value="3"/>
        <value key="Sem revestimento (leito natural)" value="4"/>
        <value key="Desconhecido" value="95"/>
        <value key="Outros" value="99"/>
        <value key="N�o" value="0"/>
        <value key="Sim" value="1"/>
        <value key="Desconhecido" value="95"/>
        <value key="Planejada" value="1"/>
        <value key="Constru�da" value="2"/>
        <value key="Abandonada" value="3"/>
        <value key="Destru�da" value="4"/>
        <value key="Em constru��o" value="5"/>
        <value key="Constru�da, mas em obras" value="6"/>
        <value key="Desconhecida" value="95"/>
        <value key="N�o aplic�vel" value="97"/>
        <value key="Peri�dico" value="2"/>
        <value key="Permanente" value="3"/>
        <value key="Tempor�rio" value="4"/>
        <value key="Desconhecido" value="95"/>
        <value key="N�o" value="0"/>
        <value key="Sim" value="1"/>
        <value key="N�o" value="0"/>
        <value key="Sim" value="1"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="tipotrechorod">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Auto-estrada" value="1"/>
        <value key="Rodovia" value="3"/>
        <value key="Entroncamento" value="4"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="jurisdicao">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Internacional" value="1"/>
        <value key="Propriedade particular" value="2"/>
        <value key="Federal" value="3"/>
        <value key="Estadual/Distrital" value="4"/>
        <value key="Municipal" value="5"/>
        <value key="Desconhecida" value="95"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="administracao">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Federal" value="1"/>
        <value key="Estadual/Distrital" value="2"/>
        <value key="Municipal" value="3"/>
        <value key="Concessionada" value="4"/>
        <value key="Privada" value="5"/>
        <value key="N�o aplic�vel" value="6"/>
        <value key="Desconhecida" value="95"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="revestimento">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Revestimento prim�rio (solto)" value="1"/>
        <value key="Pavimentado" value="2"/>
        <value key="Madeira" value="3"/>
        <value key="Sem revestimento (leito natural)" value="4"/>
        <value key="Desconhecido" value="95"/>
        <value key="Outros" value="99"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="operacional">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="N�o" value="0"/>
        <value key="Sim" value="1"/>
        <value key="Desconhecido" value="95"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="situacaofisica">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Planejada" value="1"/>
        <value key="Constru�da" value="2"/>
        <value key="Abandonada" value="3"/>
        <value key="Destru�da" value="4"/>
        <value key="Em constru��o" value="5"/>
        <value key="Constru�da, mas em obras" value="6"/>
        <value key="Desconhecida" value="95"/>
        <value key="N�o aplic�vel" value="97"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="trafego">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Peri�dico" value="2"/>
        <value key="Permanente" value="3"/>
        <value key="Tempor�rio" value="4"/>
        <value key="Desconhecido" value="95"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="trechoemperimetrourbano">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="N�o" value="0"/>
        <value key="Sim" value="1"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="tipopavimentacao">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Asfalto" value="2"/>
        <value key="Placa de concreto" value="3"/>
        <value key="Pedra regular" value="4"/>
        <value key="Ladrilho de concreto" value="5"/>
        <value key="Paralelep�pedo" value="6"/>
        <value key="Pedra irregular" value="7"/>
        <value key="Desconhecido" value="95"/>
        <value key="N�o aplic�vel" value="97"/>
        <value key="Outros" value="99"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="canteirodivisorio">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="N�o" value="0"/>
        <value key="Sim" value="1"/>
        <value key="Desconhecido" value="95"/>
      </widgetv2config>
    </edittype> 
  </edittypes>
</qgis>