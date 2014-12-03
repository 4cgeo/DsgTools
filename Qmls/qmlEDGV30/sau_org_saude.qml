<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'><qgis version="2.6.0-Brighton" minimumScale="1" maximumScale="1" simplifyDrawingHints="0" minLabelScale="0" maxLabelScale="1e+08" simplifyDrawingTol="1" simplifyMaxScale="1" hasScaleBasedVisibilityFlag="0" simplifyLocal="1" scaleBasedLabelVisibilityFlag="0"> 
  <edittypes> 
     <edittype widgetv2type="TextEdit" name="OGC_FID"> 
      <widgetv2config IsMultiline="0" fieldEditable="0" UseHtml="0" labelOnTop="0"/> 
    </edittype> 
    <edittype widgetv2type="TextEdit" name="id"> 
      <widgetv2config IsMultiline="0" fieldEditable="0" UseHtml="0" labelOnTop="0"/> 
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
    <edittype widgetv2type="ValueMap" name="grupoativecon">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Ensino m�dio" value="1"/>
        <value key="Educa��o infantil e ensino fundamental" value="2"/>
        <value key="Servi�os veterin�rios " value="3"/>
        <value key="Ensino superior" value="4"/>
        <value key="Educa��o profissional e outras atividades de ensino" value="5"/>
        <value key="Administra��o do estado e da pol�tica econ�mica e social" value="6"/>
        <value key="Servi�os coletivos prestados pela administra��o" value="7"/>
        <value key="Seguridade social" value="9"/>
        <value key="Atividades de aten��o � sa�de" value="10"/>
        <value key="Servi�o social" value="19"/>
        <value key="Desconhecido" value="95"/>
        <value key="Misto" value="98"/>
        <value key="Outros" value="99"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="classeativecon">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Atendimento hospitalar (hospital)" value="27"/>
        <value key="Atendimento �s urg�ncias e emerg�ncias (pronto-socorro)" value="28"/>
        <value key="Aten��o ambulatorial (posto e centro de sa�de)" value="29"/>
        <value key="Servi�os de complementa��o diagn�stica ou terap�utica" value="30"/>
        <value key="Outras atividades relacionadas com aten��o � sa�de (instituto de pesquisa)" value="32"/>
        <value key="Servi�os veterin�rios" value="36"/>
        <value key="Desconhecido" value="95"/>
        <value key="Misto" value="97"/>
        <value key="Outros" value="99"/>
      </widgetv2config>
    </edittype> 
  </edittypes>
</qgis>