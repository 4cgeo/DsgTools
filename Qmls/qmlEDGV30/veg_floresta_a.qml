<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'><qgis version="2.6.0-Brighton" minimumScale="1" maximumScale="1" simplifyDrawingHints="0" minLabelScale="0" maxLabelScale="1e+08" simplifyDrawingTol="1" simplifyMaxScale="1" hasScaleBasedVisibilityFlag="0" simplifyLocal="1" scaleBasedLabelVisibilityFlag="0"> 
  <edittypes> 
     <edittype widgetv2type="TextEdit" name="OGC_FID"> 
      <widgetv2config IsMultiline="0" fieldEditable="0" UseHtml="0" labelOnTop="0"/> 
    </edittype> 
    <edittype widgetv2type="TextEdit" name="id"> 
      <widgetv2config IsMultiline="0" fieldEditable="0" UseHtml="0" labelOnTop="0"/> 
    </edittype>
    <edittype widgetv2type="TextEdit" name="id_area_verde"> 
      <widgetv2config IsMultiline="0" fieldEditable="0" UseHtml="0" labelOnTop="0"/> 
    </edittype>
    <edittype widgetv2type="ValueMap" name="geometriaaproximada">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="N�o" value="0"/>
        <value key="Sim" value="1"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="tipoveg">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Floresta" value="2"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="classificacaoporte">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Arb�rea" value="4"/>
        <value key="Arbustiva" value="5"/>
        <value key="Desconhecida" value="95"/>
        <value key="Mista" value="97"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="densidade">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Alta" value="1"/>
        <value key="Baixa" value="2"/>
        <value key="Desconhecida" value="95"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="especiepredominante">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Cip�" value="11"/>
        <value key="Bambu" value="12"/>
        <value key="Arauc�ria" value="13"/>
        <value key="Sororoca" value="17"/>
        <value key="Palmeira" value="27"/>
        <value key="Sem predomin�ncia" value="37"/>
        <value key="Desconhecido" value="95"/>
        <value key="N�o identificado" value="96"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="antropizada">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="N�o" value="0"/>
        <value key="Sim" value="1"/>
        <value key="Desconhecido" value="95"/>
      </widgetv2config>
    </edittype> 
  </edittypes>
</qgis>